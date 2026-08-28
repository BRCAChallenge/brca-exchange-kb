"""
Export per-variant HerediClassify input JSON files from the database.

For each BRCA1/BRCA2 variant with GRCh38 coordinates, gathers genomic
coordinates, gnomAD v4.1 joint frequencies, BayesDel and SpliceAI scores,
and exLOVD multifactorial-likelihood data from the database, queries the
local VEP server for per-transcript consequences, and writes one JSON file
per variant conforming to HerediClassify's API/schema_input.json.

Output layout (under --out-dir):
    input/<VRS_Digest with ':' replaced by '_'>.json
    export_errors.tsv    (variants that could not be exported, with reason)

Each JSON carries an extra "_meta" key (VRS_Digest, HGVS_cDNA) ignored by
HerediClassify but used downstream by run_herediclassify.py for its summary.

Existing input files are skipped unless --overwrite is set, so reruns only
pay VEP cost for new variants.
"""

import json
import logging
import os

import click
import psycopg2
import requests

VEP_BATCH_SIZE = 200   # variants per VEP call, matches run_vep_analysis.py

log = logging.getLogger(__name__)

# BayesDel bounds from HerediClassify's schema_input.json; values outside
# would fail their input validation, so such scores are omitted.
_BAYESDEL_MIN = -1.29334
_BAYESDEL_MAX = 0.75731

_QUERY = """
SELECT v."VRS_Digest", v."Gene_Symbol", v."HGVS_cDNA",
       gc.hgvs, gc.chr, gc.pos, gc.ref, gc.alt,
       rg."Allele_frequency", rg."Allele_count",
       rg.faf95_popmax, rg.faf95_popmax_population, rg.populations,
       ab."BayesDel_nsfp33a_noAF",
       asp.result,
       ex."Combined_Prior_P", ex."Co_Occurrence_LR",
       ex."Segregation_LR", ex."Product_Of_LRs"
FROM variant v
JOIN variant_genomic_coordinates gc
     ON gc."VRS_Digest" = v."VRS_Digest" AND gc.assembly = 'GRCh38'
LEFT JOIN report_gnomad rg
     ON rg."VRS_Digest" = v."VRS_Digest"
     AND rg.version = 'v4.1' AND rg.data_type = 'joint'
LEFT JOIN analysis_bayesdel ab ON ab."VRS_Digest" = v."VRS_Digest"
LEFT JOIN analysis_spliceai asp ON asp."VRS_Digest" = v."VRS_Digest"
LEFT JOIN variant_exlovd ex ON ex."VRS_Digest" = v."VRS_Digest"
WHERE v."Gene_Symbol" = ANY(%(genes)s)
"""


def digest_to_filename(vrs_digest):
    """Filesystem-safe filename for a VRS digest (colon is legal on Linux
    but breaks many tools)."""
    return vrs_digest.replace(':', '_') + '.json'


def _defined(value):
    """Database text fields use '-' (and sometimes NULL) for 'no data'."""
    return value is not None and value != '-'


def _to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _first_int(value):
    """Parse VEP exon/intron notation ('19/23', '19-20/23') to its first number."""
    if not value:
        return None
    return _to_int(str(value).split('/')[0].split('-')[0])


def build_variant_effect(transcript_consequences, gene):
    """Map VEP transcript_consequences to HerediClassify's variant_effect list.

    Only Ensembl transcripts of the target gene with an HGVS cDNA are usable
    (schema requires transcript ~ ENST..., hgvs_c ~ c.\\S+).
    """
    effects = []
    for tc in transcript_consequences or []:
        transcript_id = tc.get('transcript_id') or ''
        if not transcript_id.startswith('ENST'):
            continue
        if gene and tc.get('gene_symbol') and tc['gene_symbol'] != gene:
            continue
        hgvsc = tc.get('hgvsc')
        if not hgvsc or ':' not in hgvsc:
            continue
        hgvs_c = hgvsc.split(':', 1)[1]
        hgvsp = tc.get('hgvsp')
        hgvs_p = hgvsp.split(':', 1)[1] if hgvsp and ':' in hgvsp else None
        consequence_terms = tc.get('consequence_terms') or []
        if not consequence_terms:
            continue
        effects.append({
            'transcript': transcript_id,
            'hgvs_c': hgvs_c,
            'hgvs_p': hgvs_p,
            'variant_type': consequence_terms,
            'exon': _first_int(tc.get('exon')),
            'intron': _first_int(tc.get('intron')),
        })
    return effects


def build_gnomad(allele_frequency, allele_count, faf95_popmax,
                 faf95_popmax_population, populations):
    """Map a report_gnomad v4.1 joint row to HerediClassify's gnomAD block.

    Returns None if the variant has no gnomAD row at all (all columns NULL
    from the LEFT JOIN) — HerediClassify treats a missing gnomAD key as
    'absent from gnomAD' (frequency 0, subpopulation "None"), which is the
    intended semantics.
    """
    if all(v is None for v in (allele_frequency, allele_count, faf95_popmax,
                               faf95_popmax_population, populations)):
        return None

    af = _to_float(allele_frequency) if _defined(allele_frequency) else None
    ac = _to_int(allele_count) if _defined(allele_count) else None

    ac_hom = 0
    popmax_af, popmax_ac = None, None
    if populations:
        for pop_data in populations.values():
            ac_hom += _to_int(pop_data.get('ac_hom')) or 0
        # popmax over the reported genetic ancestry groups
        best = max(populations.values(),
                   key=lambda p: _to_float(p.get('af')) or 0.0)
        popmax_af = _to_float(best.get('af'))
        popmax_ac = _to_int(best.get('ac'))

    faf = _to_float(faf95_popmax) if _defined(faf95_popmax) else None
    # Schema pattern for subpopulation is [A-Z]{3}; "ALL" is HerediClassify's
    # marker for "no subpopulation entry".
    subpopulation = (faf95_popmax_population.upper()
                     if _defined(faf95_popmax_population) else 'ALL')

    # All these keys are required by the schema whenever the gnomAD block is
    # present; 0 matches HerediClassify's own default for a missing key.
    return {
        'AF': af if af is not None else 0.0,
        'AC': ac if ac is not None else 0,
        'AC_hom': ac_hom,
        'subpopulation': subpopulation,
        'popmax_AF': popmax_af if popmax_af is not None else 0.0,
        'popmax_AC': popmax_ac if popmax_ac is not None else 0,
        'faf_popmax_AF': faf if faf is not None else 0.0,
    }


def build_input_json(row, vep_records):
    """Build the HerediClassify input dict for one DB row + its VEP records.

    Returns (input_dict, None) on success or (None, reason) on failure.
    """
    (vrs_digest, gene, hgvs_cdna, _hgvs, chr_, pos, ref, alt,
     allele_frequency, allele_count, faf95_popmax, faf95_popmax_population,
     populations, bayesdel, spliceai_result,
     prior, co_occurrence, segregation, product_of_lrs) = row

    if isinstance(vep_records, dict) and 'error' in vep_records:
        return None, f"VEP error: {vep_records['error']}"
    if not vep_records:
        return None, 'no VEP records'

    record = vep_records[0]
    most_severe = record.get('most_severe_consequence')
    variant_effect = build_variant_effect(
        record.get('transcript_consequences'), gene)
    if not variant_effect:
        return None, 'no usable transcript consequences'
    if not most_severe:
        most_severe = variant_effect[0]['variant_type'][0]

    data = {
        '_meta': {'VRS_Digest': vrs_digest, 'HGVS_cDNA': hgvs_cdna},
        'chr': str(chr_),
        'pos': int(pos),
        'gene': gene,
        'ref': ref,
        'alt': alt,
        'variant_type': [most_severe],
        'variant_effect': variant_effect,
    }

    gnomad = build_gnomad(allele_frequency, allele_count, faf95_popmax,
                          faf95_popmax_population, populations)
    if gnomad is not None:
        data['gnomAD'] = gnomad

    bayesdel_val = _to_float(bayesdel) if _defined(bayesdel) else None
    if bayesdel_val is not None:
        if _BAYESDEL_MIN <= bayesdel_val <= _BAYESDEL_MAX:
            data['pathogenicity_prediction_tools'] = {'BayesDel': bayesdel_val}
        else:
            log.warning('%s: BayesDel %s outside schema bounds, omitted',
                        vrs_digest, bayesdel_val)

    spliceai_val = _to_float(spliceai_result) if _defined(spliceai_result) else None
    if spliceai_val is not None and 0.0 <= spliceai_val <= 1.0:
        data['splicing_prediction_tools'] = {'SpliceAI': spliceai_val}

    for key, value in [('prior', prior),
                       ('co-occurrence', co_occurrence),
                       ('segregation', segregation),
                       ('multifactorial_log-likelihood', product_of_lrs)]:
        val = _to_float(value) if _defined(value) else None
        if val is not None:
            data[key] = val

    return data, None


def _query_vep_batch(session, vep_url, hgvs_list):
    """POST a batch of HGVS strings to the VEP server (same endpoint as
    run_vep_analysis.py). Returns dict {hgvs: [records]} or {hgvs: {'error': ...}}."""
    resp = session.post(f'{vep_url}/vep/hgvs/batch', json=hgvs_list, timeout=600)
    resp.raise_for_status()
    return resp.json()


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--vep-url', default='http://localhost:8888', show_default=True,
              envvar='VEP_SERVER_URL', help='Base URL of the local VEP REST server')
@click.option('--out-dir', required=True, type=click.Path(),
              help='Run directory; input JSONs go in <out-dir>/input/')
@click.option('--genes', default='BRCA1,BRCA2', show_default=True,
              help='Comma-separated gene symbols to export')
@click.option('--vrs-digest', default=None, metavar='DIGEST',
              help='Export only the single variant with this VRS digest')
@click.option('--limit', default=None, type=int,
              help='Export at most this many variants')
@click.option('--overwrite', is_flag=True, default=False,
              help='Rewrite input JSONs that already exist')
@click.option('--debug', is_flag=True, default=False,
              help='Print the generated JSON for each variant')
def main(db_url, schema, vep_url, out_dir, genes, vrs_digest, limit,
         overwrite, debug):
    logging.basicConfig(level=logging.WARNING, format='%(levelname)s %(message)s')

    input_dir = os.path.join(out_dir, 'input')
    os.makedirs(input_dir, exist_ok=True)
    errors_path = os.path.join(out_dir, 'export_errors.tsv')

    gene_list = [g.strip() for g in genes.split(',') if g.strip()]

    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        query, params = _QUERY, {'genes': gene_list}
        if vrs_digest:
            query += ' AND v."VRS_Digest" = %(digest)s'
            params['digest'] = vrs_digest
        with conn.cursor() as cur:
            cur.execute(query, params)
            rows = cur.fetchall()
    finally:
        conn.close()

    if vrs_digest and not rows:
        print(f'Error: VRS digest {vrs_digest!r} not found in the database.')
        return

    if not overwrite:
        rows = [r for r in rows
                if not os.path.exists(os.path.join(input_dir, digest_to_filename(r[0])))]
    if limit is not None:
        rows = rows[:limit]

    total = len(rows)
    print(f'Exporting HerediClassify input for {total} variant(s) to {input_dir} ...')

    session = requests.Session()
    written, errors = 0, []

    for batch_start in range(0, total, VEP_BATCH_SIZE):
        batch = rows[batch_start:batch_start + VEP_BATCH_SIZE]
        hgvs_map = {r[3]: r for r in batch}   # genomic hgvs -> row

        if (batch_start // VEP_BATCH_SIZE) % 10 == 0 or batch_start + VEP_BATCH_SIZE >= total:
            print(f'  {min(batch_start + VEP_BATCH_SIZE, total)}/{total}'
                  f'  ({written} written, {len(errors)} errors)')

        no_hgvs = [r for r in batch if not _defined(r[3])]
        for r in no_hgvs:
            errors.append((r[0], 'no GRCh38 genomic HGVS'))

        hgvs_list = [h for h in hgvs_map if _defined(h)]
        if not hgvs_list:
            continue
        try:
            results = _query_vep_batch(session, vep_url, hgvs_list)
        except Exception as e:
            log.warning('VEP batch starting at %d failed: %s', batch_start, e)
            errors.extend((hgvs_map[h][0], f'VEP batch failed: {e}') for h in hgvs_list)
            continue

        for hgvs, records in results.items():
            row = hgvs_map.get(hgvs)
            if row is None:
                continue
            data, reason = build_input_json(row, records)
            if data is None:
                errors.append((row[0], reason))
                continue
            if debug:
                print(json.dumps(data, indent=2))
            with open(os.path.join(input_dir, digest_to_filename(row[0])), 'w') as f:
                json.dump(data, f, indent=2)
            written += 1

    with open(errors_path, 'a') as f:
        for digest, reason in errors:
            f.write(f'{digest}\t{reason}\n')

    print(f'Done. Wrote {written} input file(s), {len(errors)} error(s)'
          f'{" (see " + errors_path + ")" if errors else ""}.')


if __name__ == '__main__':
    main()

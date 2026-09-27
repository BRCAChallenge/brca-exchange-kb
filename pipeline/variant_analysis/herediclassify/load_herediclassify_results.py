"""
Load a HerediClassify run's per-variant JSON results into analysis_herediclassify.

Reads <run-dir>/output/*.json as written by variant_analysis/run_herediclassify.py
and upserts one row per variant. Safe to re-run while a run is still in progress:
the upsert keys on (VRS_Digest, method_name), so a later pass picks up newly
finished variants and refreshes any that were reclassified.

Usage:
    python load_herediclassify_results.py \
        --run-dir /data/new_schema/data_out/herediclassify_runs/2026-09-07 \
        --method-name herediclassify_1.1 --schema pipeline

Everything describing how a result was produced goes into provenance_metadata:
the run label, the gnomAD release the frequency rules were fed, the tool and
config versions, and the classification timestamp. method_name is a generated
column derived from that JSON, so it is never written directly.

Runs before 2026-09-07 predate the per-variant gnomAD provenance fields. That
export fed gnomAD v4.1 joint to every variant, so --assume-gnomad-release lets
that be recorded rather than lost; newer runs carry it per variant and need no
flag.
"""

import argparse
import glob
import json
import os
import sys

import psycopg2
import psycopg2.extras

DB_URL = 'postgresql://postgres:postgres@localhost/storage.pg'
BATCH = 1000

# Copied from the result JSON into provenance_metadata when present.
_PROVENANCE_KEYS = ('popfreq_code', 'gnomad_version', 'gnomad_data_type',
                    'config_name', 'config_version',
                    'herediclassify_version', 'herediclassify_commit')

_UPSERT = """
INSERT INTO analysis_herediclassify
    ("VRS_Digest", provenance_metadata, classification_protein,
     classification_splicing, rules)
VALUES %s
ON CONFLICT ("VRS_Digest", method_name) DO UPDATE
    SET provenance_metadata     = EXCLUDED.provenance_metadata,
        classification_protein  = EXCLUDED.classification_protein,
        classification_splicing = EXCLUDED.classification_splicing,
        rules                   = EXCLUDED.rules
"""


def build_row(result, method_name, assume_gnomad):
    """Turn one result JSON into the tuple the upsert expects.

    Returns None if the file carries no VRS digest to key on.
    """
    digest = result.get('VRS_Digest')
    if not digest:
        return None

    provenance = {'method_name': method_name}
    for key in _PROVENANCE_KEYS:
        if result.get(key) is not None:
            provenance[key] = result[key]
    if assume_gnomad and 'gnomad_version' not in provenance:
        version, data_type = assume_gnomad
        provenance['gnomad_version'] = version
        provenance['gnomad_data_type'] = data_type
    if result.get('timestamp'):
        provenance['classified_at'] = result['timestamp']

    rules = result.get('rules') or {}
    # classification_protein/_splicing sit alongside the rules in the tool's
    # output; lift them into their own columns and leave the rest as the rules
    # object, which holds only the per-rule dicts.
    rule_dicts = {k: v for k, v in rules.items() if isinstance(v, dict)}

    return (digest,
            json.dumps(provenance),
            rules.get('classification_protein'),
            rules.get('classification_splicing'),
            json.dumps(rule_dicts))


def load(run_dir, db_url, schema, method_name, assume_gnomad, dry_run, limit):
    paths = sorted(glob.glob(os.path.join(run_dir, 'output', '*.json')))
    if limit:
        paths = paths[:limit]
    print(f'Reading {len(paths)} result file(s) from {run_dir}/output ...')

    rows, malformed = [], 0
    for path in paths:
        try:
            with open(path) as fh:
                result = json.load(fh)
        except (ValueError, OSError) as e:
            malformed += 1
            print(f'  skipping {os.path.basename(path)}: {e}', file=sys.stderr)
            continue
        row = build_row(result, method_name, assume_gnomad)
        if row is None:
            malformed += 1
            continue
        rows.append(row)

    print(f'Prepared {len(rows)} row(s)' + (f', {malformed} unreadable' if malformed else ''))
    if dry_run:
        if rows:
            digest, provenance, prot, splice, _ = rows[0]
            print('Dry run; first row would be:')
            print(f'  VRS_Digest={digest}')
            print(f'  provenance_metadata={provenance}')
            print(f'  classification_protein={prot} classification_splicing={splice}')
        return

    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        # Digests absent from variant would violate the foreign key and abort
        # the whole load, so drop them up front and say how many.
        with conn.cursor() as cur:
            cur.execute('SELECT "VRS_Digest" FROM variant')
            known = {r[0] for r in cur}
        unknown = [r for r in rows if r[0] not in known]
        if unknown:
            print(f'Skipping {len(unknown)} row(s) whose VRS digest is not in variant')
            rows = [r for r in rows if r[0] in known]

        with conn.cursor() as cur:
            for i in range(0, len(rows), BATCH):
                psycopg2.extras.execute_values(cur, _UPSERT, rows[i:i + BATCH])
                print(f'  {min(i + BATCH, len(rows))}/{len(rows)} written')
        conn.commit()
        with conn.cursor() as cur:
            cur.execute('SELECT count(*) FROM analysis_herediclassify WHERE method_name = %s',
                        (method_name,))
            total = cur.fetchone()[0]
        print(f'Done. {schema}.analysis_herediclassify now holds {total} row(s) '
              f'for method {method_name!r}.')
    finally:
        conn.close()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run-dir', required=True,
                   help='HerediClassify run directory containing output/')
    p.add_argument('--db-url', default=os.environ.get('PIPELINE_DB_URL', DB_URL))
    p.add_argument('--schema', default='pipeline')
    p.add_argument('--method-name', default='herediclassify_1.0',
                   help='Run label; stored in provenance_metadata and surfaced '
                        'as the generated method_name column')
    p.add_argument('--assume-gnomad-release', default=None, metavar='VERSION:TYPE',
                   help='Record this gnomAD release for results that predate the '
                        'per-variant provenance fields, e.g. v4.1:joint')
    p.add_argument('--limit', type=int, default=None)
    p.add_argument('--dry-run', action='store_true')
    args = p.parse_args()

    assume_gnomad = None
    if args.assume_gnomad_release:
        if ':' not in args.assume_gnomad_release:
            p.error('--assume-gnomad-release must look like v4.1:joint')
        assume_gnomad = tuple(args.assume_gnomad_release.split(':', 1))

    load(args.run_dir, args.db_url, args.schema, args.method_name,
         assume_gnomad, args.dry_run, args.limit)


if __name__ == '__main__':
    main()

"""
Load an ARIANE run's per-variant JSON results into analysis_ariane.

Reads <out-dir>/results/*.json as written by variant_analysis/query_ariane.py
and upserts one row per variant. Safe to re-run: a later pass picks up variants
classified since, and refreshes any that were reclassified.

Usage:
    python load_ariane_results.py \
        --out-dir /data/new_schema/data_out/ariane \
        --method-name ariane_1.9 --schema vcep_development

provenance_metadata carries the run label plus what ARIANE reported about how
the result was produced: its application and build versions (surfaced as the
generated ariane_version/build_version columns), the policy applied, the API
version, the classification engine and fingerprint, and the generation time.
The response itself is stored whole -- result.criteria in criteria, everything
else in result_detail.
"""

import argparse
import collections
import glob
import json
import os
import sys

import psycopg2
import psycopg2.extras

DB_URL = 'postgresql://postgres:postgres@localhost/storage.pg'
BATCH = 1000

# provenance_metadata key -> key in _meta.ariane as ARIANE reports it.
_ARIANE_PROVENANCE = {
    'ariane_version':           'application_version',
    'build_version':            'build_version',
    'build_version_checked_at': 'build_version_checked_at',
    'api_version':              'api_version',
    'classification_engine':    'classification_engine',
    'classifier_fingerprint':   'classifier_fingerprint',
    'request_id':               'request_id',
    'classified_at':            'generated_at',
}

_UPSERT = """
INSERT INTO analysis_ariane
    ("VRS_Digest", provenance_metadata, criteria, result_detail)
VALUES %s
ON CONFLICT ("VRS_Digest", method_name) DO UPDATE
    SET provenance_metadata = EXCLUDED.provenance_metadata,
        criteria            = EXCLUDED.criteria,
        result_detail       = EXCLUDED.result_detail
"""


def build_row(payload, method_name):
    """Turn one result JSON into the tuple the upsert expects.

    Returns None if the file carries no VRS digest or no classification.
    """
    meta = payload.get('_meta') or {}
    result = payload.get('result')
    digest = meta.get('VRS_Digest')
    if not digest or not result:
        return None

    ariane = meta.get('ariane') or {}
    provenance = {'method_name': method_name}
    for key, source_key in _ARIANE_PROVENANCE.items():
        if ariane.get(source_key) is not None:
            provenance[key] = ariane[source_key]
    if ariane.get('policy') is not None:
        provenance['policy'] = ariane['policy']
    for key in ('source', 'retrieved', 'gene', 'c_notation'):
        if meta.get(key) is not None:
            provenance[key] = meta[key]

    criteria = result.get('criteria') or []
    detail = {k: v for k, v in result.items() if k != 'criteria'}
    return (digest, json.dumps(provenance), json.dumps(criteria), json.dumps(detail))


def load(out_dir, db_url, schema, method_name, dry_run, limit):
    paths = sorted(glob.glob(os.path.join(out_dir, 'results', '*.json')))
    if limit:
        paths = paths[:limit]
    print(f'Reading {len(paths)} result file(s) from {out_dir}/results ...')

    rows, malformed = [], 0
    versions, builds, labels = collections.Counter(), collections.Counter(), collections.Counter()
    for path in paths:
        try:
            with open(path) as fh:
                payload = json.load(fh)
        except (ValueError, OSError) as e:
            malformed += 1
            print(f'  skipping {os.path.basename(path)}: {e}', file=sys.stderr)
            continue
        row = build_row(payload, method_name)
        if row is None:
            malformed += 1
            continue
        rows.append(row)
        provenance = json.loads(row[1])
        versions[provenance.get('ariane_version')] += 1
        builds[provenance.get('build_version')] += 1
        labels[(payload['result'] or {}).get('predicted_label')] += 1

    print(f'Prepared {len(rows)} row(s)' + (f', {malformed} unreadable' if malformed else ''))
    print(f'  ariane_version: {dict(versions)}')
    print(f'  build_version:  {dict(builds)}')
    print(f'  predicted_label: {dict(labels)}')
    if dry_run:
        if rows:
            digest, provenance, criteria, detail = rows[0]
            print('Dry run; first row would be:')
            print(f'  VRS_Digest={digest}')
            print(f'  provenance_metadata={provenance}')
            print(f'  criteria: {len(json.loads(criteria))} entries, '
                  f'result_detail: {len(json.loads(detail))} keys')
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
            cur.execute('SELECT count(*) FROM analysis_ariane WHERE method_name = %s',
                        (method_name,))
            total = cur.fetchone()[0]
        print(f'Done. {schema}.analysis_ariane now holds {total} row(s) '
              f'for method {method_name!r}.')
    finally:
        conn.close()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out-dir', required=True,
                   help='ARIANE run directory containing results/')
    p.add_argument('--db-url', default=os.environ.get('PIPELINE_DB_URL', DB_URL))
    p.add_argument('--schema', default='pipeline')
    p.add_argument('--method-name', default='ariane_1.9',
                   help='Run label; stored in provenance_metadata and surfaced '
                        'as the generated method_name column')
    p.add_argument('--limit', type=int, default=None)
    p.add_argument('--dry-run', action='store_true')
    args = p.parse_args()

    load(args.out_dir, args.db_url, args.schema, args.method_name,
         args.dry_run, args.limit)


if __name__ == '__main__':
    main()

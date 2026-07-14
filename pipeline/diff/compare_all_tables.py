#!/usr/bin/env python
"""
Compare gnomAD-related tables (variant_gnomad and report_gnomad) across two schemas.

Usage:
    python compare_all_tables.py --old OLD_SCHEMA --new NEW_SCHEMA --output-dir DIR

Runs compare_table.py for each table and writes output files to DIR:
    variant_gnomad.detail.txt
    report_gnomad.detail.txt

report_gnomad uses (VRS_Digest_id, version, data_type) as the natural key for
matching rows, since its surrogate id column differs across schemas.
"""

import argparse
import os
import subprocess
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_COMPARE = os.path.join(_SCRIPT_DIR, 'compare_table.py')

TABLES = [
    dict(table='variant_gnomad', pk=None),
    dict(table='report_gnomad',  pk=['VRS_Digest_id', 'version', 'data_type']),
]


def run(old, new, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    python = sys.executable
    for spec in TABLES:
        table = spec['table']
        output = os.path.join(output_dir, f'{table}.detail.txt')
        cmd = [python, _COMPARE,
               '--table', table,
               '--old', old,
               '--new', new,
               '--output', output]
        if spec['pk']:
            cmd += ['--pk'] + spec['pk']
        print(f'\n=== {table} ===')
        subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--old',        required=True, help='Old schema name')
    parser.add_argument('--new',        required=True, help='New schema name')
    parser.add_argument('--output-dir', required=True, help='Directory for output files')
    args = parser.parse_args()
    run(args.old, args.new, args.output_dir)


if __name__ == '__main__':
    main()

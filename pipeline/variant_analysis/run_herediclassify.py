#!/usr/bin/env python3
"""
Run HerediClassify over exported input JSONs, writing results to the filesystem.

IMPORTANT: this script must run under HerediClassify's own virtualenv (it
imports HerediClassify's classify() in-process, which needs pyensembl,
cyvcf2, etc.). It deliberately imports nothing from this pipeline — stdlib
plus HerediClassify's own venv packages — so the venv split between stage 1
(export_herediclassify_input.py, pipeline venv) and stage 2 (this script)
cannot bite.

Reads <run-dir>/input/*.json (produced by export_herediclassify_input.py),
calls HerediClassify's classify() for each, and writes:

    <run-dir>/output/<same filename>.json   rule results + final classes + metadata
    <run-dir>/classify_errors.tsv           variants whose classification raised
    <run-dir>/herediclassify_summary.tsv    one row per classified variant

Existing outputs are skipped unless --overwrite, so an interrupted run is
resumable. Per-variant exceptions are logged and do not stop the run.
"""

import argparse
import csv
import datetime
import glob
import json
import os
import pathlib
import subprocess
import sys


def get_commit(herediclassify_dir):
    try:
        return subprocess.run(
            ['git', '-C', herediclassify_dir, 'rev-parse', 'HEAD'],
            capture_output=True, text=True, timeout=30,
        ).stdout.strip() or None
    except Exception:
        return None


def load_classify(herediclassify_dir):
    """Import HerediClassify's classify() (and its version) in-process.

    HerediClassify's modules use flat sibling imports, so its
    variant_classification directory must be on sys.path. Importing
    classify also loads pyensembl (Ensembl release 110) — slow once,
    then reused for every variant.
    """
    # pybedtools needs `bedtools` on PATH; the venv's bin dir carries a
    # symlink to it, so make sure the running interpreter's bin dir is first.
    os.environ['PATH'] = (os.path.dirname(sys.executable) + os.pathsep
                          + os.environ.get('PATH', ''))
    sys.path.insert(0, os.path.join(herediclassify_dir, 'variant_classification'))
    from classify import classify          # noqa: E402
    import pybedtools                      # noqa: E402
    try:
        from _version import __version__   # noqa: E402
    except ImportError:
        __version__ = None
    return classify, __version__, pybedtools


def summarise(output_dir, summary_path):
    """Write one TSV row per classified variant from the output JSONs."""
    outputs = []
    for path in sorted(glob.glob(os.path.join(output_dir, '*.json'))):
        with open(path) as f:
            outputs.append(json.load(f))
    if not outputs:
        print('No outputs to summarise.')
        return

    rule_names = sorted({name
                         for out in outputs
                         for name in out['rules']
                         if not name.startswith('classification_')})

    with open(summary_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['VRS_Digest', 'gene', 'HGVS_cDNA',
                         'classification_protein', 'classification_splicing']
                        + rule_names)
        for out in outputs:
            rules = out['rules']
            row = [out.get('VRS_Digest'), out.get('gene'), out.get('HGVS_cDNA'),
                   rules.get('classification_protein'),
                   rules.get('classification_splicing')]
            for name in rule_names:
                rule = rules.get(name)
                if rule is None:
                    row.append('')
                elif rule.get('status'):
                    row.append(f"met({rule.get('strength')})")
                else:
                    row.append('not_met')
            writer.writerow(row)
    print(f'Summary: {len(outputs)} variant(s) -> {summary_path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--herediclassify-dir', required=True,
                        help='Path to the HerediClassify git clone')
    parser.add_argument('--config', default=None,
                        help='HerediClassify config yaml '
                             '(default: <herediclassify-dir>/config.yaml)')
    parser.add_argument('--run-dir', required=True,
                        help='Run directory containing input/ from '
                             'export_herediclassify_input.py')
    parser.add_argument('--vrs-digest', default=None,
                        help='Classify only this variant (matched against the '
                             'input JSON _meta.VRS_Digest or filename)')
    parser.add_argument('--limit', type=int, default=None,
                        help='Classify at most this many variants')
    parser.add_argument('--overwrite', action='store_true',
                        help='Re-classify variants that already have an output JSON')
    parser.add_argument('--summary-only', action='store_true',
                        help='Skip classification, just rebuild the summary TSV '
                             'from existing outputs')
    args = parser.parse_args()

    config_path = args.config or os.path.join(args.herediclassify_dir, 'config.yaml')
    input_dir = os.path.join(args.run_dir, 'input')
    output_dir = os.path.join(args.run_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    errors_path = os.path.join(args.run_dir, 'classify_errors.tsv')
    summary_path = os.path.join(args.run_dir, 'herediclassify_summary.tsv')

    if args.summary_only:
        summarise(output_dir, summary_path)
        return

    input_paths = sorted(glob.glob(os.path.join(input_dir, '*.json')))
    if args.vrs_digest:
        sanitized = args.vrs_digest.replace(':', '_')
        input_paths = [p for p in input_paths
                       if os.path.basename(p) == sanitized + '.json']
        if not input_paths:
            sys.exit(f'Error: no input JSON for VRS digest {args.vrs_digest!r} '
                     f'in {input_dir}')
    if not args.overwrite:
        input_paths = [p for p in input_paths
                       if not os.path.exists(
                           os.path.join(output_dir, os.path.basename(p)))]
    if args.limit is not None:
        input_paths = input_paths[:args.limit]

    total = len(input_paths)
    print(f'Classifying {total} variant(s) with config {config_path} ...')
    if total == 0:
        summarise(output_dir, summary_path)
        return

    print('Loading HerediClassify (pyensembl etc., this takes a while) ...')
    classify, hc_version, pybedtools = load_classify(args.herediclassify_dir)
    commit = get_commit(args.herediclassify_dir)

    # classify() leaks pybedtools temp files (freed only atexit), which can
    # fill the filesystem over a long run — keep them off /tmp and delete
    # them per variant, as HerediClassify's own webservice does.
    pybedtools_tmp = os.path.join(args.run_dir, 'pybedtools_tmp')
    os.makedirs(pybedtools_tmp, exist_ok=True)
    pybedtools.set_tempdir(pybedtools_tmp)

    done, failed = 0, 0
    started = datetime.datetime.now()
    with open(errors_path, 'a') as errors_fh:
        for i, path in enumerate(input_paths, start=1):
            with open(path) as f:
                data = json.load(f)
            meta = data.get('_meta', {})
            digest = meta.get('VRS_Digest') or os.path.basename(path)[:-len('.json')]

            try:
                final_config, result_str = classify(
                    pathlib.Path(config_path), json.dumps(data))
                rules = json.loads(result_str)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                failed += 1
                errors_fh.write(f'{digest}\t{type(e).__name__}\t{e}\n')
                errors_fh.flush()
                continue
            finally:
                pybedtools.cleanup()

            out = {
                'VRS_Digest': digest,
                'gene': data.get('gene'),
                'HGVS_cDNA': meta.get('HGVS_cDNA'),
                'config_name': final_config.get('name'),
                'config_version': final_config.get('version'),
                'herediclassify_version': hc_version,
                'herediclassify_commit': commit,
                'timestamp': datetime.datetime.now().isoformat(timespec='seconds'),
                'rules': rules,
            }
            with open(os.path.join(output_dir, os.path.basename(path)), 'w') as f:
                json.dump(out, f, indent=2)
            done += 1

            if i % 100 == 0 or i == total:
                elapsed = (datetime.datetime.now() - started).total_seconds()
                rate = elapsed / i if i else 0
                print(f'  {i}/{total}  ({done} ok, {failed} failed, '
                      f'{rate:.2f}s/variant)')

    print(f'Done. {done} classified, {failed} failed'
          f'{" (see " + errors_path + ")" if failed else ""}.')
    summarise(output_dir, summary_path)


if __name__ == '__main__':
    main()

# HerediClassify installation and run guide

[HerediClassify](https://github.com/akatzke/HerediClassify) (GPL-3.0) is run
**whole** as an external tool in its own clone + virtualenv — no code from it
is vendored into this repository (a few one-line patches were applied to the
clone itself; see "Local patches" below). The one file copied here is its
`API/schema_input.json` (kept beside this doc, provenance noted inside), so
the export test can validate against it without a clone. The bridge is two
scripts with a filesystem handoff:

| Stage | Script | Venv | Does |
|---|---|---|---|
| 1 | `pipeline/variant_analysis/export_herediclassify_input.py` | `~mcline/.venv/brcaexchange` | DB + VEP server → one input JSON per variant |
| 2 | `pipeline/variant_analysis/run_herediclassify.py` | HerediClassify's venv | input JSONs → classify() → output JSONs + summary TSV |

Nothing is written to the database. Output root convention:
`/data/new_schema/data_out/herediclassify_runs/<YYYY-MM-DD>/`.

## As-built installation record (2026-08-24)

Install root: `/data/new_schema/herediclassify_tool/`
Clone: `/data/new_schema/herediclassify_tool/HerediClassify`
commit `22905fad9602b31582e89960435cc5611818cede`, HerediClassify version 1.0.6.

### 1. Clone + venv

Built with **Python 3.12.3** (the box has no 3.10; their pins install fine on
3.12 with the tweaks below):

```bash
cd /data/new_schema/herediclassify_tool
git clone https://github.com/akatzke/HerediClassify
cd HerediClassify
python3.12 -m venv venv
venv/bin/pip install --upgrade pip
# requirements.patched.txt = requirements.txt minus dev-only pins
# (nose, pylint, ipython toolchain, mock, openpyxl, yoyo-migrations ...)
# nose in particular does not install on modern Python.
venv/bin/pip install -r requirements.patched.txt   # (pybedtools removed too)
# pybedtools 0.9.1 needs both of these to build:
venv/bin/pip install 'setuptools<81'        # pkg_resources still needed (hgvs too)
venv/bin/pip install --no-build-isolation pybedtools==0.9.1
```

### 2. bedtools

No system bedtools/samtools and no sudo; static binary instead:

```bash
mkdir -p /data/new_schema/herediclassify_tool/bin
curl -sL -o bin/bedtools \
  https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
chmod +x bin/bedtools
ln -s /data/new_schema/herediclassify_tool/bin/bedtools \
      HerediClassify/venv/bin/bedtools
```

`run_herediclassify.py` prepends its interpreter's bin dir to PATH so
pybedtools finds this symlink. samtools was not needed.

### 3. Reference data

- **pyensembl / Ensembl 110**: already cached at `~/.cache/pyensembl/GRCh38/ensembl110`
  (~1 GB, GTF db prebuilt) — no download was needed; the FASTA `.pickle`
  indexes were built on first use. If it ever disappears:
  `venv/bin/pyensembl install --release 110 --species human`.
- **ClinVar** (185 MB): `databases/Clinvar/clinvar.vcf.gz` from
  `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/`, then
  `venv/bin/python install_dependencies/data_filter_clinvar.py -i .../clinvar.vcf.gz -f true`
  → `clinvar_snv_two_star.vcf.gz`, `clinvar_small_indel_two_star.vcf.gz`
  (+ three_star variants), all tabix-indexed.
- **Uniprot repeats**: UCSC `unipRepeat` track JSON →
  `install_dependencies/data_format_uniprot_rep.py` →
  `databases/Uniprot/repeats_hg38_uniprot.bed`.
- **MANE v1.3**: the `current/` FTP dir no longer carries v1.3 — use
  `https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz`,
  then `install_dependencies/data_format_MANE_transcript_list.py` →
  `databases/MANE/MANE.GRCh38.v1.3.ensembl_genomic.gtf_transcript_list.csv`.
- The critical-region BEDs, splice-site tables, and PM5 annotations ship in
  the clone's `data/` directory (no download).

### 4. Config edits (in the clone)

- `config.yaml` and `gene_specific/acmg_brca1.yaml` / `acmg_brca2.yaml`:
  `annotation_files.root: "/data/new_schema/herediclassify_tool/"`;
  `config.yaml` `gene_specific_configs.root:
  "/data/new_schema/herediclassify_tool/HerediClassify/gene_specific/"`.
  (Gene-specific configs fully *replace* the top-level config at runtime.)
- **Dropped `ps1_protein_enigma` and `ps1_splicing_clingen`** from both BRCA
  gene configs. These need `clinvar_snv_spliceai` — a ClinVar VCF annotated
  with SpliceAI scores that upstream's own installer no longer builds (their
  `merge_clinvar_spliceai*.sh` needs ngs-bits plus Illumina's
  auth-gated precomputed SpliceAI VCFs). With a plain VCF the annotation
  step raises (missing `SpliceAI` column in `format_spliceai`), killing every
  variant, so dropping the rules is the only safe degradation. PS1 is
  therefore **never met** in this run configuration — revisit if a
  SpliceAI-annotated ClinVar (even BRCA-region-scoped) is built later.

### 5. Local patches to the clone (upstream bugs)

- `variant_classification/acmg_rules/__init__.py`: commented out
  `from acmg_rules.pp4_pcd import *` — the module is referenced but not
  committed upstream.
- `variant_classification/config_annotation.py:114`: commented out the
  `"pp4_pcd": Rules.Pp4_pcd` mapping for the same reason.

### 6. Smoke test (passed 2026-08-24)

```bash
cd /data/new_schema/herediclassify_tool/HerediClassify
PATH=/data/new_schema/herediclassify_tool/bin:$PATH venv/bin/python \
    variant_classification/classify.py -c config.yaml \
    -i "$(cat test/test_variants_gene_specific/BRCA1_frameshift_variant.json | tr -d '\n')"
# → full rule dict, classification_protein: 5 (Pathogenic)
```

## Running

```bash
RUN_DIR=/data/new_schema/data_out/herediclassify_runs/$(date +%F)

# Stage 1 — pipeline venv (VEP server must be up, see pipeline/vep_server)
~mcline/.venv/brcaexchange/bin/python \
    pipeline/variant_analysis/export_herediclassify_input.py --out-dir $RUN_DIR

# Stage 2 — HerediClassify venv (~0.13 s/variant after a slow first load)
/data/new_schema/herediclassify_tool/HerediClassify/venv/bin/python \
    pipeline/variant_analysis/run_herediclassify.py \
    --herediclassify-dir /data/new_schema/herediclassify_tool/HerediClassify \
    --run-dir $RUN_DIR
```

Both stages skip already-written files, so both are resumable / incremental;
use `--overwrite` to force. `--vrs-digest ga4gh:VA...` runs a single variant;
stage 2 `--summary-only` rebuilds the TSV from existing outputs; `--limit N`
caps a run.

Or via Luigi (`ExportHerediClassifyInput` / `RunHerediClassify` in
`pipeline/workflow/variant_analysis.py`), which defaults the run dir to the
convention above.

## Output layout

```
<run-dir>/
  input/<VRS_Digest, ':'→'_'>.json    stage 1: HerediClassify input JSON
  output/<same name>.json             stage 2: rules + final classes + metadata
  export.log                          stage 1 console log (when run detached)
  export_errors.tsv                   stage 1 failures (digest, reason)
  classify_errors.tsv                 stage 2 failures (digest, exception, message)
  herediclassify_summary.tsv          per-variant: gene, HGVS, final protein/splicing
                                      class (ACMG 1-5), met/not-met per rule
```

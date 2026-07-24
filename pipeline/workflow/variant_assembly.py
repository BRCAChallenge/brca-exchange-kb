import csv
import datetime
import os
import re
import shutil
import tempfile

import luigi
from luigi.util import requires

luigi.auto_namespace(scope=__name__)

from common import config as brca_config
from workflow import pipeline_utils
from workflow.pipeline_common import PipelineParams

_pipeline_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_RESOURCES_DIR = os.path.normpath(os.path.join(_pipeline_dir, '..', '..', 'resources'))
clinvar_method_dir           = os.path.join(_pipeline_dir, 'clinvar')
functional_assays_method_dir = os.path.join(_pipeline_dir, 'functional_assays')
data_merging_method_dir      = os.path.join(_pipeline_dir, 'data_merging')

# GRCh38 chromosome lengths — keyed by both bare (13) and prefixed (chr13) forms
_GRCh38_CONTIGS = {
    '1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
    '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
    '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
    '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
    '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
    '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415,
    'MT': 16569,
}
_GRCh38_CONTIGS.update({f'chr{k}': v for k, v in _GRCh38_CONTIGS.items()})


class VCFAssemblyTask(luigi.Task):
    """Base task for variant_assembly. Creates only the directories this pipeline uses."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cfg = PipelineParams.get_instance()

        file_parent_dir = os.path.abspath(self.cfg.file_parent_dir)
        output_dir      = os.path.abspath(self.cfg.output_dir)
        create = pipeline_utils.create_path_if_nonexistent
        self.artifacts_dir    = create(output_dir      + "/release/artifacts")
        self.clinvar_file_dir = create(file_parent_dir + "/ClinVar")
        self.ex_lovd_file_dir = create(file_parent_dir + "/exLOVD")
        self.lovd_file_dir    = create(file_parent_dir + "/LOVD")
        self.enigma_file_dir  = create(file_parent_dir + "/enigma")
        self.assays_dir       = create(file_parent_dir + "/functional_assays")
        self.gnomad_file_dir  = create(file_parent_dir + "/gnomAD")
        self.vcf_dir          = create(file_parent_dir + "/vcf")
        self.varaico_dir      = create(file_parent_dir + "/varaico")

    def on_failure(self, exception):
        def _rename_file(path):
            if os.path.exists(path):
                ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                failed_name = f"FAILED_{ts}_{os.path.basename(path)}"
                os.rename(path, os.path.join(os.path.dirname(path), failed_name))

        outputs = self.output()
        if isinstance(outputs, luigi.LocalTarget):
            _rename_file(outputs.path)
        else:
            for o in outputs.values():
                _rename_file(o.path)
        return super().on_failure(exception)

    def _normalize_vcf_header(self, input_vcf, output_vcf):
        """Copy input_vcf to output_vcf, injecting missing ##contig and ##INFO lines."""
        header, data = [], []
        with open(input_vcf) as fh:
            for line in fh:
                if line.startswith('#'):
                    header.append(line)
                else:
                    data.append(line)

        # Inject missing ##contig lines
        declared_contigs = set()
        for line in header:
            m = re.search(r'^##contig=<ID=([^,>]+)', line)
            if m:
                declared_contigs.add(m.group(1))

        referenced_contigs = {row.split('\t')[0] for row in data if row.strip()}
        missing_contigs = [c for c in (referenced_contigs - declared_contigs)
                           if c in _GRCh38_CONTIGS]

        # Inject missing ##INFO lines
        declared_info = set()
        for line in header:
            m = re.search(r'^##INFO=<ID=([^,>]+)', line)
            if m:
                declared_info.add(m.group(1))

        referenced_info = set()
        for row in data:
            if not row.strip():
                continue
            cols = row.rstrip('\n').split('\t')
            if len(cols) > 7 and cols[7] not in ('.', ''):
                for field in cols[7].split(';'):
                    key = field.split('=')[0]
                    if key:
                        referenced_info.add(key)
        missing_info = sorted(referenced_info - declared_info)

        if not missing_contigs and not missing_info:
            with open(output_vcf, 'w') as fh:
                fh.writelines(header)
                fh.writelines(data)
            return

        chrom_idx = next(i for i, l in enumerate(header) if l.startswith('#CHROM'))
        inject = []
        if missing_contigs:
            inject += [f'##contig=<ID={c},length={_GRCh38_CONTIGS[c]},assembly=GRCh38>\n'
                       for c in sorted(missing_contigs)]
        if missing_info:
            inject += [f'##INFO=<ID={k},Number=.,Type=String,Description="{k}">\n'
                       for k in missing_info]
        header = header[:chrom_idx] + inject + header[chrom_idx:]
        with open(output_vcf, 'w') as fh:
            fh.writelines(header)
            fh.writelines(data)

    def _run_process_with_pipeline_path(self, args, **kwargs):
        """Run a subprocess with _pipeline_dir prepended to PYTHONPATH."""
        prev = os.environ.get('PYTHONPATH', '')
        os.environ['PYTHONPATH'] = _pipeline_dir + ((':' + prev) if prev else '')
        try:
            pipeline_utils.run_process(args, **kwargs)
        finally:
            os.environ['PYTHONPATH'] = prev

    def _vrs_annotate(self, input_vcf, output_vcf_gz, output_pkl):
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as tmp:
            tmp_path = tmp.name
        try:
            self._normalize_vcf_header(input_vcf, tmp_path)
            args = ["vrs-annotate", "vcf", tmp_path,
                    "--vcf-out", output_vcf_gz,
                    "--pkl-out", output_pkl,
                    "--dataproxy-uri", f"seqrepo+file://{self.cfg.seq_repo_dir}"]
            pipeline_utils.run_process(args)
            pipeline_utils.check_file_for_contents(output_vcf_gz)
        finally:
            os.unlink(tmp_path)


###############################################
#                   CLINVAR                   #
###############################################

class DownloadLatestClinvarData(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVarVCVRelease_00-latest.xml.gz")

    def run(self):
        os.chdir(self.clinvar_file_dir)
        clinvar_data_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarVCVRelease_00-latest.xml.gz"
        pipeline_utils.download_file_and_display_progress(clinvar_data_url)


@requires(DownloadLatestClinvarData)
class ConvertLatestClinvarDataToXML(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVar.xml")

    def run(self):
        os.chdir(clinvar_method_dir)
        genes_opts = [s for g in self.cfg.gene_metadata['symbol'] for s in ['--gene', g]]
        pipeline_utils.run_process(["python", "filter_clinvar.py",
                                    "-i", self.input().path,
                                    "-o", self.output().path] + genes_opts)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertLatestClinvarDataToXML)
class ConvertClinvarXMLToTXT(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVar.txt")

    def run(self):
        os.chdir(clinvar_method_dir)
        args = ["python", "clinVarParse.py",
                self.input().path,
                "--logs", self.clinvar_file_dir + "/clinvar_xml_to_txt.log",
                "--assembly", "GRCh38"]
        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertClinvarXMLToTXT)
class ConvertClinvarTXTToVCF(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(self.clinvar_file_dir + "/ClinVar.vcf")

    def run(self):
        os.chdir(data_merging_method_dir)
        args = ["python", "convert_tsv_to_vcf.py",
                "-i", self.input().path,
                "-o", self.output().path,
                "-s", "ClinVar"]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertClinvarTXTToVCF)
class VRSAnnotateClinVar(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/ClinVar.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/ClinVar.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


###############################################
#                  exLOVD                     #
###############################################

class DownloadStaticExLOVDVCF(VCFAssemblyTask):
    exlovd_static_vcf_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/exLOVD.BRCA12.sorted.hg38.vcf',
        description='URL to download static exLOVD vcf from')

    def output(self):
        return luigi.LocalTarget(f"{self.ex_lovd_file_dir}/exLOVD.BRCA12.sorted.hg38.vcf")

    def run(self):
        data = pipeline_utils.urlopen_with_retry(self.exlovd_static_vcf_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)


@requires(DownloadStaticExLOVDVCF)
class VRSAnnotateEXLOVD(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/exLOVD.BRCA12.sorted.hg38.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/exLOVD.BRCA12.sorted.hg38.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


@requires(VRSAnnotateEXLOVD)
class DeduplicateExLOVD(VCFAssemblyTask):
    """Remove duplicate exLOVD records that share the same VRS digest.

    Keeps the record citing the most recent year in key_observational_reference;
    ties broken by number of non-missing INFO fields.
    The pkl is unchanged (already keyed by VRS digest, so inherently unique).
    """

    def output(self):
        annotated = self.input()
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/exLOVD.BRCA12.sorted.hg38.dedup.vcf.gz"),
            "pkl": annotated["pkl"],
        }

    def run(self):
        script = os.path.join(_pipeline_dir, "variant_assembly", "deduplicate_exlovd_vcf.py")
        args = [
            "python", script,
            "--input-vcf",  self.input()["vcf"].path,
            "--output-vcf", self.output()["vcf"].path,
        ]
        self._run_process_with_pipeline_path(args)


##############################################
#               sharedLOVD                   #
##############################################

class DownloadStaticSharedLOVDVCF(VCFAssemblyTask):
    shared_lovd_static_vcf_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/LOVD.sorted.hg38.vcf',
        description='URL to download static sharedLOVD vcf from')

    def output(self):
        return luigi.LocalTarget(f"{self.lovd_file_dir}/LOVD.sorted.hg38.vcf")

    def run(self):
        data = pipeline_utils.urlopen_with_retry(self.shared_lovd_static_vcf_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)


@requires(DownloadStaticSharedLOVDVCF)
class VRSAnnotateSharedLOVD(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/LOVD.sorted.hg38.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/LOVD.sorted.hg38.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


###############################################
#                  ENIGMA                     #
###############################################

@requires(ConvertLatestClinvarDataToXML)
class FilterEnigmaAssertions(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.enigma_file_dir, 'enigma_clinvar.xml'))

    def run(self):
        os.chdir(clinvar_method_dir)
        args = ["python", "filter_enigma_data.py",
                self.input().path,
                self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(FilterEnigmaAssertions)
class ExtractEnigmaFromClinvar(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.enigma_file_dir, 'enigma.tsv'))

    def run(self):
        os.chdir(clinvar_method_dir)
        genes_opts = [s for g in self.cfg.gene_metadata['symbol'] for s in ['--gene', g]]
        args = ["python", "enigma_from_clinvar.py",
                self.input().path,
                self.output().path,
                '--logs', os.path.join(self.enigma_file_dir, 'enigma.log'),
               ] + genes_opts
        self._run_process_with_pipeline_path(args)


@requires(ExtractEnigmaFromClinvar)
class ConvertEnigmaToVCF(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.enigma_file_dir, 'enigma.vcf'))

    def run(self):
        os.chdir(data_merging_method_dir)
        args = ["python", "convert_tsv_to_vcf.py",
                "-i", self.input().path,
                "-o", self.output().path,
                "-s", "ENIGMA"]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertEnigmaToVCF)
class VRSAnnotateEnigma(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/enigma.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/enigma.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


###############################################
#             FUNCTIONAL ASSAYS               #
###############################################

class DownloadFunctionalAssaysInputFile(VCFAssemblyTask):
    findlay_BRCA1_ring_function_scores_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/ENIGMA_BRCA12_FunctionalAssays_latest_with_function_scores.tsv',
        description='URL to download functional assay data from')

    def output(self):
        return luigi.LocalTarget(
            self.assays_dir + "/ENIGMA_BRCA12_FunctionalAssays_latest_with_function_scores.tsv")

    def run(self):
        data = pipeline_utils.urlopen_with_retry(
            self.findlay_BRCA1_ring_function_scores_url).read()
        with open(self.output().path, "wb") as f:
            f.write(data)


@requires(DownloadFunctionalAssaysInputFile)
class ConvertFunctionalAssaysToVCF(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(
            self.assays_dir + "/ENIGMA_BRCA12_functional_assays_scores.hg19.vcf")

    def run(self):
        os.chdir(functional_assays_method_dir)
        args = ["python", "convert_functional_assay_tsv_to_vcf.py",
                "-v", "-i", self.input().path,
                "-o", self.output().path,
                "-a", "functionalAssayAnnotation",
                "-l", self.artifacts_dir + "/functional_assays_error_variants.log",
                "-s", "ENIGMABRCA12FunctionalAssaysFunctionScores"]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ConvertFunctionalAssaysToVCF)
class CrossmapFunctionalAssays(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(
            os.path.join(self.assays_dir, "ENIGMA_BRCA12_functional_assays_scores.hg38.vcf"))

    def run(self):
        brca_resources_dir = self.cfg.resources_dir
        args = ["CrossMap.py", "vcf",
                brca_resources_dir + "/hg19ToHg38.over.chain.gz",
                self.input().path,
                brca_resources_dir + "/hg38.fa",
                self.output().path]
        pipeline_utils.run_process(args)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(CrossmapFunctionalAssays)
class SortFunctionalAssays(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(
            os.path.join(self.assays_dir, "ENIGMA_BRCA12_functional_assays_scores.sorted.hg38.vcf"))

    def run(self):
        args = ["vcf-sort", self.input().path]
        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(SortFunctionalAssays)
class VRSAnnotateFunctionalAssays(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.vcf_dir}/ENIGMA_BRCA12_functional_assays_scores.sorted.hg38.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.vcf_dir}/ENIGMA_BRCA12_functional_assays_scores.sorted.hg38.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


###############################################
#                  gnomAD                     #
###############################################

class GnomADTask(VCFAssemblyTask):
    def _download_file(self, url, output_path):
        data = pipeline_utils.urlopen_with_retry(url).read()
        with open(output_path, "wb") as f:
            f.write(data)


class DownloadStaticGnomADVCF(GnomADTask):
    gnomAD_v2_static_vcf_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/gnomADv2.sorted.hg38.vcf',
        description='URL to download static gnomAD v2 vcf from')
    gnomAD_v3_static_data_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/gnomADv3.sorted.hg38.vcf',
        description='URL to download static gnomAD v3 vcf from')
    gnomad_v4_static_data_url = luigi.Parameter(
        default='https://brcaexchange.org/backend/downloads/vcf/gnomADv4.sorted.hg38.vcf',
        description='URL to download static gnomAD v4 vcf from')

    def output(self):
        return {
            "v2": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv2.sorted.hg38.vcf"),
            "v3": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv3.sorted.hg38.vcf"),
            "v4": luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv4.sorted.hg38.vcf"),
        }

    def run(self):
        self._download_file(self.gnomAD_v2_static_vcf_url, self.output()["v2"].path)
        self._download_file(self.gnomAD_v3_static_data_url, self.output()["v3"].path)
        self._download_file(self.gnomad_v4_static_data_url, self.output()["v4"].path)


@requires(DownloadStaticGnomADVCF)
class PostprocessGnomADv3VCF(VCFAssemblyTask):
    """Add coverage INFO field to the downloaded gnomAD v3 VCF.

    The static v3 VCF already contains variant_id and flags fields but lacks
    coverage. This task adds coverage from the v3.1 genome parquet so that
    LoadVCFsToDatabase can store per-variant coverage for v3 entries.
    """

    gnomad_v3_genome_coverage = luigi.Parameter(
        default=os.path.join(_RESOURCES_DIR, 'gnomADv3.1.coverage.genome.parquet'),
        description='Path to the gnomAD v3.1 genome coverage parquet')

    def output(self):
        return luigi.LocalTarget(
            f"{self.gnomad_file_dir}/gnomADv3.sorted.hg38.postprocessed.vcf"
        )

    def run(self):
        download_out = self.input()
        script = os.path.join(_pipeline_dir, 'gnomad', 'download_gnomad_fourpointone.py')
        args = [
            'python', script,
            '--postprocess-only',
            '--input-vcf', download_out['v3'].path,
            '-c', self.gnomad_v3_genome_coverage,
            '-o', self.output().path,
            '-v',
        ]
        self._run_process_with_pipeline_path(args)


@requires(DownloadStaticGnomADVCF, PostprocessGnomADv3VCF)
class VRSAnnotateGnomAD(VCFAssemblyTask):
    def output(self):
        return {
            "v2_vcf": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv2.sorted.hg38.vcf.gz"),
            "v2_pkl": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv2.sorted.hg38.dicts.pkl"),
            "v3_vcf": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv3.sorted.hg38.vcf.gz"),
            "v3_pkl": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv3.sorted.hg38.dicts.pkl"),
            "v4_vcf": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv4.sorted.hg38.vcf.gz"),
            "v4_pkl": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv4.sorted.hg38.dicts.pkl"),
        }

    def run(self):
        download_out, v3_postprocessed = self.input()
        self._vrs_annotate(
            download_out["v2"].path,
            self.output()["v2_vcf"].path,
            self.output()["v2_pkl"].path,
        )
        self._vrs_annotate(
            v3_postprocessed.path,
            self.output()["v3_vcf"].path,
            self.output()["v3_pkl"].path,
        )
        self._vrs_annotate(
            download_out["v4"].path,
            self.output()["v4_vcf"].path,
            self.output()["v4_pkl"].path,
        )


class DownloadGnomADCoverage(VCFAssemblyTask):
    """Download gnomAD genome + exome coverage and write a combined weighted-mean parquet.

    Also writes per-source parquets for the V4 exome and V3 genome coverage.
    Disabled by default; pass --DownloadGnomADCoverage-enabled true to run.
    """
    enabled = luigi.BoolParameter(default=False)
    gene_config = luigi.Parameter(
        default=os.path.join(_pipeline_dir, 'workflow', 'gene_config_brca_only.txt'),
        description='Gene config file with chr/start_hg38/end_hg38 columns')
    coverage_output = luigi.Parameter(
        default=os.path.join(_pipeline_dir, '..', '..', 'resources', 'gnomADv4.1.coverage.joint.parquet'),
        description='Output path for the combined weighted-mean coverage parquet')
    exome_coverage_output = luigi.Parameter(
        default=os.path.join(_pipeline_dir, '..', '..', 'resources', 'gnomADv4.1.coverage.exome.parquet'),
        description='Output path for the V4 exome gene-region coverage parquet')
    genome_coverage_output = luigi.Parameter(
        default=os.path.join(_pipeline_dir, '..', '..', 'resources', 'gnomADv3.1.coverage.genome.parquet'),
        description='Output path for the V3 genome gene-region coverage parquet')

    def output(self):
        return {
            'joint':  luigi.LocalTarget(os.path.normpath(self.coverage_output)),
            'exome':  luigi.LocalTarget(os.path.normpath(self.exome_coverage_output)),
            'genome': luigi.LocalTarget(os.path.normpath(self.genome_coverage_output)),
        }

    def run(self):
        if not self.enabled:
            return
        script = os.path.join(_pipeline_dir, 'gnomad', 'download_coverage_v4.py')
        args = [
            'python', script,
            '-g', self.gene_config,
            '-c', os.path.normpath(self.coverage_output),
            '--exome-coverage-output', os.path.normpath(self.exome_coverage_output),
            '--genome-coverage-output', os.path.normpath(self.genome_coverage_output),
            '-v',
        ]
        self._run_process_with_pipeline_path(args)


@requires(DownloadGnomADCoverage)
class ProvideGnomADv41JointVCF(GnomADTask):
    """Provide the gnomAD v4.1 joint VCF for BRCA genes.

    Copies from the resources directory if available; otherwise downloads
    and postprocesses the data from gnomAD, then copies the result back to
    the resources directory for future runs.
    """
    resources_vcf = luigi.Parameter(
        default=os.path.join(_RESOURCES_DIR, 'gnomADv4.1.joint.hg38.vcf'),
        description='Cached VCF path in the resources directory')

    def output(self):
        return luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv4.1.joint.hg38.vcf")

    def run(self):
        try:
            shutil.copy(self.resources_vcf, self.output().path)
        except (FileNotFoundError, OSError):
            script = os.path.join(_pipeline_dir, 'gnomad', 'download_gnomad_fourpointone.py')
            args = [
                'python', script,
                '-g', self.gene_config,
                '-c', self.input()['joint'].path,
                '-o', self.output().path,
                '--source', 'joint',
                '-v',
            ]
            self._run_process_with_pipeline_path(args)
            shutil.copy(self.output().path, self.resources_vcf)


@requires(DownloadGnomADCoverage)
class ProvideGnomADv41ExomeVCF(GnomADTask):
    """Provide the gnomAD v4.1 exome VCF for BRCA genes.

    Copies from the resources directory if available; otherwise downloads
    and postprocesses the data from gnomAD, then copies the result back to
    the resources directory for future runs.
    """
    resources_vcf = luigi.Parameter(
        default=os.path.join(_RESOURCES_DIR, 'gnomADv4.1.exome.hg38.vcf'),
        description='Cached VCF path in the resources directory')

    def output(self):
        return luigi.LocalTarget(f"{self.gnomad_file_dir}/gnomADv4.1.exome.hg38.vcf")

    def run(self):
        try:
            shutil.copy(self.resources_vcf, self.output().path)
        except (FileNotFoundError, OSError):
            script = os.path.join(_pipeline_dir, 'gnomad', 'download_gnomad_fourpointone.py')
            args = [
                'python', script,
                '-g', self.gene_config,
                '-c', self.input()['exome'].path,
                '-o', self.output().path,
                '--source', 'exome',
                '-v',
            ]
            self._run_process_with_pipeline_path(args)
            shutil.copy(self.output().path, self.resources_vcf)


@requires(ProvideGnomADv41JointVCF, ProvideGnomADv41ExomeVCF)
class VRSAnnotateGnomADv41(VCFAssemblyTask):
    """VRS-annotate the gnomAD v4.1 joint and exome VCFs."""

    def output(self):
        return {
            "joint_vcf": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv4.1.joint.hg38.vcf.gz"),
            "joint_pkl": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv4.1.joint.hg38.dicts.pkl"),
            "exome_vcf": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv4.1.exome.hg38.vcf.gz"),
            "exome_pkl": luigi.LocalTarget(f"{self.vcf_dir}/gnomADv4.1.exome.hg38.dicts.pkl"),
        }

    def run(self):
        joint_in, exome_in = self.input()
        self._vrs_annotate(joint_in.path,
                           self.output()["joint_vcf"].path,
                           self.output()["joint_pkl"].path)
        # gnomAD v4.1 exome ships with VRS_Allele_IDs already in the header;
        # strip them so vrs-annotate can add its own.
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as tmp:
            stripped = tmp.name
        try:
            pipeline_utils.run_process([
                'bcftools', 'annotate', '-x',
                'INFO/VRS_Allele_IDs,INFO/VRS_Starts,INFO/VRS_Ends,INFO/VRS_States',
                '-o', stripped, exome_in.path
            ])
            self._vrs_annotate(stripped,
                               self.output()["exome_vcf"].path,
                               self.output()["exome_pkl"].path)
        finally:
            os.unlink(stripped)


###############################################
#           LOAD VCFs TO DATABASE             #
###############################################

@requires(
    VRSAnnotateEnigma,
    VRSAnnotateClinVar,
    VRSAnnotateSharedLOVD,
    DeduplicateExLOVD,
    VRSAnnotateGnomAD,
    VRSAnnotateFunctionalAssays,
    VRSAnnotateGnomADv41,
)
class LoadVCFsToDatabase(VCFAssemblyTask):
    """Call the load_vcf Django management command to load all VCFs into the pipeline DB."""

    django_dir = luigi.Parameter(
        default='/data/new_schema/code/website/django',
        description='Directory containing Django manage.py')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, "load_vcfs_to_db.done"))

    def run(self):
        enigma_in, clinvar_in, lovd_in, exlovd_in, gnomad_in, assays_in, gnomad_v41_in = self.input()
        gene_config = os.path.join(_pipeline_dir, 'workflow', 'gene_config_brca_only.txt')
        args = [
            "python", "manage.py", "load_vcf", "--skip-checks",
            "--gene-config",             gene_config,
            "--enigma-vcf",              enigma_in["vcf"].path,
            "--enigma-pkl",              enigma_in["pkl"].path,
            "--clinvar-vcf",             clinvar_in["vcf"].path,
            "--clinvar-pkl",             clinvar_in["pkl"].path,
            "--lovd-vcf",                lovd_in["vcf"].path,
            "--lovd-pkl",                lovd_in["pkl"].path,
            "--exlovd-vcf",              exlovd_in["vcf"].path,
            "--exlovd-pkl",              exlovd_in["pkl"].path,
            "--gnomad-v2-vcf",           gnomad_in["v2_vcf"].path,
            "--gnomad-v2-pkl",           gnomad_in["v2_pkl"].path,
            "--gnomad-v3-vcf",           gnomad_in["v3_vcf"].path,
            "--gnomad-v3-pkl",           gnomad_in["v3_pkl"].path,
            "--gnomad-v41-joint-vcf",    gnomad_v41_in["joint_vcf"].path,
            "--gnomad-v41-joint-pkl",    gnomad_v41_in["joint_pkl"].path,
            "--gnomad-v41-exome-vcf",    gnomad_v41_in["exome_vcf"].path,
            "--gnomad-v41-exome-pkl",    gnomad_v41_in["exome_pkl"].path,
            "--functional-assay-vcf",    assays_in["vcf"].path,
            "--functional-assay-pkl",    assays_in["pkl"].path,
        ]
        os.chdir(self.django_dir)
        os.environ['PIPELINE_SCHEMA'] = self.cfg.db_schema
        pipeline_utils.run_process(args)
        with open(self.output().path, "w") as f:
            f.write("done\n")


###############################################
#         VARAICO LITERATURE MINING           #
###############################################

# NOTE: requires `bigBedToBed` in PATH, same as `bedToBigBed` in genomeBrowserTrack/ — download
# from https://hgdownload.soe.ucsc.edu/admin/exe/ . It supports remote random-access queries via
# -chrom/-start/-end against a URL, so the ~475MB varaico.bb is never downloaded in full.
_varaico_bb_url = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/varaico.bb"


class ExtractVaraicoBRCARegions(VCFAssemblyTask):
    """Extract BRCA1/BRCA2 rows from the UCSC varaico bigBed track via remote random access."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.varaico_dir, "varaico_brca_raw.tsv"))

    def run(self):
        gene_config = brca_config.load_config(
            os.path.join(_pipeline_dir, 'workflow', 'gene_config_brca_only.txt'))

        header_written = False
        with open(self.output().path, "w") as out_f:
            for _, gene in gene_config.iterrows():
                chrom = f"chr{gene['chr']}"
                tmp_path = self.output().path + f".{chrom}.tmp"
                args = [
                    "bigBedToBed",
                    f"-chrom={chrom}",
                    f"-start={int(gene['start_hg38'])}",
                    f"-end={int(gene['end_hg38'])}",
                    "-tsv",
                    _varaico_bb_url,
                    tmp_path,
                ]
                pipeline_utils.run_process(args)
                with open(tmp_path) as tmp_f:
                    lines = tmp_f.readlines()
                os.remove(tmp_path)
                if not lines:
                    continue
                if not header_written:
                    out_f.write(lines[0])
                    header_written = True
                out_f.writelines(lines[1:])
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(ExtractVaraicoBRCARegions)
class FilterVaraicoToGeneBoundaries(VCFAssemblyTask):
    """Re-validate gene-boundary containment and build a minimal VCF of the distinct variants."""

    def output(self):
        return {
            "tsv": luigi.LocalTarget(os.path.join(self.varaico_dir, "varaico_mentions_filtered.tsv")),
            "vcf": luigi.LocalTarget(os.path.join(self.varaico_dir, "varaico_variants.vcf")),
        }

    def run(self):
        gene_config = brca_config.load_config(
            os.path.join(_pipeline_dir, 'workflow', 'gene_config_brca_only.txt'))
        boundaries = {
            f"chr{row['chr']}": (int(row['start_hg38']), int(row['end_hg38']))
            for _, row in gene_config.iterrows()
        }

        with open(self.input().path, newline='') as in_f:
            reader = csv.DictReader(in_f, delimiter='\t', quoting=csv.QUOTE_NONE)
            fieldnames = reader.fieldnames
            kept_rows = []
            for row in reader:
                bounds = boundaries.get(row['chrom'])
                if bounds is None:
                    continue
                gene_start, gene_end = bounds
                if int(row['chromStart']) >= gene_start and int(row['chromEnd']) <= gene_end:
                    kept_rows.append(row)

        with open(self.output()["tsv"].path, "w", newline='') as out_f:
            writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(kept_rows)

        # Distinct (chrom, pos, ref, alt) combos only — VRS digests don't depend on which
        # paper/snippet a row came from, and the same variant is frequently mentioned by
        # multiple papers. ID = "chrom:pos:ref:alt" so the DB loader can recompute the same
        # key from each raw row without a separate lookup file.
        distinct_variants = {}
        for row in kept_rows:
            pos = int(row['chromStart']) + 1
            key = f"{row['chrom']}:{pos}:{row['ref']}:{row['alt']}"
            distinct_variants[key] = (row['chrom'], pos, row['ref'], row['alt'])

        varaico_contigs = sorted({v[0] for v in distinct_variants.values()})
        with open(self.output()["vcf"].path, "w") as vcf_f:
            vcf_f.write("##fileformat=VCFv4.2\n")
            for contig in varaico_contigs:
                vcf_f.write(f"##contig=<ID={contig},length={_GRCh38_CONTIGS[contig]}>\n")
            vcf_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for key, (chrom, pos, ref, alt) in sorted(
                    distinct_variants.items(), key=lambda kv: (kv[1][0], kv[1][1])):
                vcf_f.write(f"{chrom}\t{pos}\t{key}\t{ref}\t{alt}\t.\t.\t.\n")


@requires(FilterVaraicoToGeneBoundaries)
class SortVaraicoVCF(VCFAssemblyTask):
    def output(self):
        return luigi.LocalTarget(os.path.join(self.varaico_dir, "varaico_variants.sorted.vcf"))

    def run(self):
        args = ["vcf-sort", self.input()["vcf"].path]
        pipeline_utils.run_process(args, redirect_stdout_path=self.output().path)
        pipeline_utils.check_file_for_contents(self.output().path)


@requires(SortVaraicoVCF)
class VRSAnnotateVaraico(VCFAssemblyTask):
    def output(self):
        return {
            "vcf": luigi.LocalTarget(f"{self.varaico_dir}/varaico_variants.sorted.vcf.gz"),
            "pkl": luigi.LocalTarget(f"{self.varaico_dir}/varaico_variants.sorted.dicts.pkl"),
        }

    def run(self):
        self._vrs_annotate(self.input().path,
                           self.output()["vcf"].path,
                           self.output()["pkl"].path)


###############################################
#         COMPUTE GENOMIC HGVS STRINGS        #
###############################################

@requires(LoadVCFsToDatabase)
class ComputeGenomicHGVS(VCFAssemblyTask):
    """Populate genomic_coordinates.hgvs with normalized genomic HGVS strings."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, "compute_genomic_hgvs.done"))

    def run(self):
        script = os.path.join(_pipeline_dir, "variant_analysis", "compute_genomic_hgvs.py")
        args = ["python", script, "--schema", self.cfg.db_schema]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, "w") as f:
            f.write("done\n")


###############################################
#        QUERY CLINGEN ALLELE REGISTRY        #
###############################################

@requires(ComputeGenomicHGVS)
class QueryClinGenAlleleRegistry(VCFAssemblyTask):
    """Populate variant CA_ID, title, HGVS_cDNA, ensembl_cdna, HGVS_Protein,
    ensembl_protein, and Synonyms from the ClinGen Allele Registry, and insert
    GRCh37 rows into genomic_coordinates."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, "query_clingen_allele_registry.done"))

    def run(self):
        script = os.path.join(_pipeline_dir, "variant_assembly", "query_clingen_allele_registry.py")
        args = ["python", script, "--schema", self.cfg.db_schema]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, "w") as f:
            f.write("done\n")


###############################################
#          RUN PSEUDONYM GENERATOR            #
###############################################

_data_merging_dir = os.path.join(_pipeline_dir, 'data_merging')

@requires(LoadVCFsToDatabase)
class RunPseudonymGenerator(VCFAssemblyTask):
    """Compute HGVS pseudonyms for all variants in the DB via pseudonym_generator.py.

    Reads variants directly from the pipeline DB and writes enriched fields
    (Reference_Sequence, HGVS_cDNA, HGVS_Protein, CA_ID, Synonyms, GRCh37
    genomic_coordinates) back to the DB."""

    db_url = luigi.Parameter(
        default='postgresql://postgres:postgres@localhost/storage.pg',
        description='PostgreSQL connection URL for the pipeline DB')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, "run_pseudonym_generator.done"))

    def run(self):
        script = os.path.join(_pipeline_dir, "variant_assembly", "pseudonym_generator.py")
        gene_config = os.path.join(_pipeline_dir, 'workflow', 'gene_config_brca_only.txt')
        args = [
            "python", script,
            "--db-url",    self.db_url,
            "--schema",    self.cfg.db_schema,
            "--configfile", gene_config,
            "--resources",  self.cfg.resources_dir,
        ]
        os.environ['HGVS_SEQREPO_DIR'] = self.cfg.seq_repo_dir
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, "w") as f:
            f.write("done\n")


###############################################
#           LOAD ENIGMA DOMAINS               #
###############################################

class LoadEnigmaDomains(VCFAssemblyTask):
    """Populate enigma_domain with ENIGMA CI functional domain coordinates."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, "load_enigma_domains.done"))

    def run(self):
        script = os.path.join(_pipeline_dir, "variant_analysis", "load_enigma_domains.py")
        args = ["python", script, "--schema", self.cfg.db_schema]
        self._run_process_with_pipeline_path(args)
        with open(self.output().path, "w") as f:
            f.write("done\n")


###############################################
#      LOAD VARAICO PAPERS (FINAL STEP)       #
###############################################

@requires(VRSAnnotateVaraico, FilterVaraicoToGeneBoundaries, LoadVCFsToDatabase,
          QueryClinGenAlleleRegistry, LoadEnigmaDomains, RunPseudonymGenerator)
class LoadVaraicoPapersToDatabase(VCFAssemblyTask):
    """Call the add_varaico_papers Django management command to populate paper /
    variant_in_paper for varaico-mined BRCA1/BRCA2 variants already present in `variant`.

    Ordered as the final step of VCFAssembly: depends on every other task that
    populates or enriches `variant`/`genomic_coordinates` (LoadVCFsToDatabase,
    QueryClinGenAlleleRegistry, LoadEnigmaDomains, RunPseudonymGenerator), not just
    LoadVCFsToDatabase, so it always runs last instead of racing them."""

    django_dir = luigi.Parameter(
        default='/data/new_schema/code/website/django',
        description='Directory containing Django manage.py')

    def output(self):
        return luigi.LocalTarget(os.path.join(self.varaico_dir, "load_varaico_papers_to_db.done"))

    def run(self):
        vrs_in, filter_in, _, _, _, _ = self.input()
        args = [
            "python", "manage.py", "add_varaico_papers", "--skip-checks",
            "--raw-tsv",           filter_in["tsv"].path,
            "--vrs-annotated-vcf", vrs_in["vcf"].path,
        ]
        os.chdir(self.django_dir)
        os.environ['PIPELINE_SCHEMA'] = self.cfg.db_schema
        pipeline_utils.run_process(args)
        with open(self.output().path, "w") as f:
            f.write("done\n")


###############################################
#               TOP-LEVEL TASK                #
###############################################

@requires(
    VRSAnnotateClinVar,
    VRSAnnotateEXLOVD,
    VRSAnnotateSharedLOVD,
    VRSAnnotateEnigma,
    VRSAnnotateFunctionalAssays,
    VRSAnnotateGnomAD,
    QueryClinGenAlleleRegistry,
    LoadEnigmaDomains,
    RunPseudonymGenerator,
    LoadVaraicoPapersToDatabase,
)
class VCFAssembly(VCFAssemblyTask):
    """Runs all source chains, VRS-annotates each VCF, loads to DB, and writes sentinel."""

    def output(self):
        return luigi.LocalTarget(os.path.join(self.vcf_dir, "variant_assembly.done"))

    def run(self):
        with open(self.output().path, "w") as f:
            f.write("done\n")

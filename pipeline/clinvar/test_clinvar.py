from . import clinvar_common

from common import hgvs_utils, variant_utils
import xml.etree.ElementTree as ET


def test_simple_genomic_coordinate_extraction():
    sample_location = """
      <Location>
        <CytogeneticLocation>13q13.1</CytogeneticLocation>
        <SequenceLocation Assembly="GRCh38" AssemblyAccessionVersion="GCF_000001405.38" forDisplay="true" AssemblyStatus="cur\
rent" Chr="13" Accession="NC_000013.11" start="32345245" stop="32345245" display_start="32345245" display_stop="32345245" var\
iantLength="1" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        <SequenceLocation Assembly="GRCh37" AssemblyAccessionVersion="GCF_000001405.25" AssemblyStatus="previous" Chr="13" Ac\
cession="NC_000013.10" start="32919382" stop="32919382" display_start="32919382" display_stop="32919382" variantLength="1" po\
sitionVCF="32919382" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
      </Location>
    """

    location_el = ET.fromstring(sample_location)

    genomic_coords = clinvar_common.extract_genomic_coordinates_from_location(location_el)

    assert genomic_coords[hgvs_utils.HgvsWrapper.GRCh38_Assem] == variant_utils.VCFVariant(13, 32345245, "A", "G")


def test_genomic_coordinate_extraction_from_NM():

    # 3'UTR position (NM_007294.3's CDS is 5592bp, so a plain c.6591 position,
    # as opposed to c.*999, is out of bounds and would raise an HGVS error).
    sample_name="NM_007294.3:c.*999_*1000del"
    genomic_coords = clinvar_common.accession_to_genomic_coordinates(sample_name)

    assert genomic_coords[hgvs_utils.HgvsWrapper.GRCh38_Assem] == variant_utils.VCFVariant(17, 43044677, "AAT", "A")


def test_preprocess_element_value():
    assert clinvar_common._preprocess_element_value('NM_000059.3(BRCA2):c.6591_6592del (p.Glu2198fs)') == 'NM_000059.3(BRCA2):c.6591_6592del'


def test_is_bic_designation():
    assert clinvar_common.is_bic_designation('1294del41')
    assert clinvar_common.is_bic_designation('5277A>G')
    assert clinvar_common.is_bic_designation('999insA')
    assert not clinvar_common.is_bic_designation('p.(Leu392GlnfsTer6)')
    assert clinvar_common.is_bic_designation('S76X')
    assert clinvar_common.is_bic_designation('R245X')
    assert clinvar_common.is_bic_designation('E1703V')
    assert clinvar_common.is_bic_designation('S309=')
    assert not clinvar_common.is_bic_designation('NM_000059.4(BRCA2):c.67G>A')
    assert not clinvar_common.is_bic_designation('NM_007294.4:c.3541del')
    assert not clinvar_common.is_bic_designation('NC_000013.10:g.32900419G>T')
    assert not clinvar_common.is_bic_designation('LRG_292t1:c.2066_2069delGTAA')


def test_gene_chromosome_mismatch_rejected():
    # GeneList tags this variant as BRCA1, but its own coordinates are on
    # chromosome 13 (BRCA2's chromosome) -- a real ClinVar data-quality issue.
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="13" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        </Location>
        <GeneList>
          <Gene Symbol="BRCA1"/>
        </GeneList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, gene_chromosomes={'BRCA1': 17, 'BRCA2': 13}, debug=False)

    assert v.geneSymbol is None


def test_gene_chromosome_match_accepted():
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="17" positionVCF="43123988" referenceAlleleVCF="C" alternateAlleleVCF="T"/>
        </Location>
        <GeneList>
          <Gene Symbol="BRCA1"/>
        </GeneList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, gene_chromosomes={'BRCA1': 17, 'BRCA2': 13}, debug=False)

    assert v.geneSymbol == 'BRCA1'


def test_gene_chromosome_check_skipped_when_not_provided():
    # Without gene_chromosomes, behavior is unchanged: no cross-check happens.
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="13" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        </Location>
        <GeneList>
          <Gene Symbol="BRCA1"/>
        </GeneList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, debug=False)

    assert v.geneSymbol == 'BRCA1'


def test_variant_bic_nomenclature_from_other_name_list():
    # OtherNameList mixes a legacy BIC-style designator with an HGVS-style
    # protein change synonym; only the former should be picked out.
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="13" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        </Location>
        <OtherNameList>
          <Name>1294del41</Name>
          <Name>p.(Leu392GlnfsTer6)</Name>
        </OtherNameList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, debug=False)

    assert v.bic_nomenclature == '1294del41'
    assert 'p.(Leu392GlnfsTer6)' in v.synonyms


def test_variant_bic_nomenclature_absent():
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="13" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        </Location>
        <OtherNameList>
          <Name>p.(Leu392GlnfsTer6)</Name>
        </OtherNameList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, debug=False)

    assert v.bic_nomenclature is None








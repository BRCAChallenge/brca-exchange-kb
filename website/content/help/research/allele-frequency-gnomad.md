
gnomAD provides _BRCA1_ and _BRCA2_ allele frequencies for the BRCA Exchange. The goal of gnomAD is to “aggregate and harmonize exome and genome sequencing data from a variety of large-scale sequencing projects, and to make summary data available for the wider scientific community.” Its predecessor, the Exome Aggregation Consortium (ExAC) contained only exome data. The gnomAD data set used by BRCA Exchange is the “non-cancer” subset, and excludes data from TCGA and other cancer cohorts. This ensures that frequencies used to represent normal, unaffected populations are not skewed by the inclusion of cancer patients. 

As an aggregation database, gnomAD data subsumes much of ExAC and 1000 Genomes data. 

#### Exomes and Genomes 
gnomAD provides frequencies from both genome (Version 3.1) and joint (Version 4.1) sequencing projects. On the BRCA Exchange, these frequencies are presented in different subtiles. A small "E" indicates Exomes data, a small "G" icon indicates Genomes data, and a small "J" icon indicates both Exomes and Genomes data. If data is available on both Exomes and Genomes, a total of four subtiles will be available:

* gnomAD V3.1 Genomes, Non-Cancer (Graphical)
* gnomAD V3.1 Genomes, Non-Cancer (Numerical)
* gnomAD V4.1 Joint (Graphical)
* gnomAD V4.1 Joint (Numerical)

This data is updated during the monthly release cycle, and will contain any new information provided by gnomAD. 

#### gnomAD Populations

The BRCA Exchange displays the primary subsets used by ExAC and gnomAD (i.e. “East  Asian”). Additional detailed population data, such as female and male distributions, can be found in gnomAD. All data are from gnomAD's non-cancer subset, meaning that these allele frequencies do not include observations from any cancer studies.

#### The following fields are available under gnomAD 3.1 Genomes, Non-Cancer (Numerical):

* #### Total Allele frequency ((Allele_frequency_genome_GnomADv3))
	Overall allele frequency (for genomes) according to all studied populations, per gnomAD 3.1
* #### Popmax Filtering Allele Frequency (95% confidence) ((faf95_popmax_genome_GnomADv3))
	The 95% confidence bound on the most significant gnomAD genome population frequency (gnomAD 3.1)
* #### Allele frequency African ((Allele_frequency_genome_AFR_GnomADv3))
	Allele Frequency in genomes of African/African American Populations (gnomAD 3.1)
* #### Allele frequency Latino ((Allele_frequency_genome_AMR_GnomADv3))
	Allele Frequency in genomes of Latino Populations (gnomAD 3.1)
* #### Allele frequency Ashkenazi Jewish ((Allele_frequency_genome_ASJ_GnomADv3))
	Allele Frequency in genomes of Ashkenazi Jewish Populations (gnomAD 3.1)
* #### Allele frequency East Asian ((Allele_frequency_genome_EAS_GnomADv3))
	Allele Frequency in genomes of East Asian Populations (gnomAD 3.1)
* #### Allele frequency Finnish ((Allele_frequency_genome_FIN_GnomADv3))
	Allele Frequency in genomes of Finnish Populations (gnomAD 3.1)
* #### Allele frequency Non-Finnish European ((Allele_frequency_genome_NFE_GnomADv3))
	Allele Frequency in genomes of Non-Finnish European Populations (gnomAD 3.1)
* #### Allele frequency Other ((Allele_frequency_genome_OTH_GnomADv3))
	Allele Frequency in genomes of Other Populations (gnomAD 3.1)
* #### Allele frequency South Asian ((Allele_frequency_genome_SAS_GnomADv3))
	Allele Frequency in genomes of South Asian Populations (gnomAD 3.1)
* #### Allele frequency Middle Eastern ((Allele_frequency_genome_MID_GnomADv3))
	Allele Frequency in genomes of Middle Eastern Populations (gnomAD 3.1)
* #### Allele frequency Amish ((Allele_frequency_genome_AMI_GnomADv3))
	Allele Frequency in genomes of Amish Populations (gnomAD 3.1)

#### The following fields are available under gnomAD v4.1 Joint (Numerical):

* #### Total Allele frequency ((Allele_frequency_joint_GnomADv41))
Overall allele frequency according to all studied populations, combining gnomAD v4.1 exome and genome data (gnomAD v4.1 Joint)
* #### Grpmax Filtering Allele Frequency (95% confidence) ((faf95_popmax_joint_GnomADv41))
The 95% confidence bound on the highest filtering allele frequency across gnomAD v4.1 genetic ancestry groups. The Amish, Ashkenazi Jewish, Finnish, Middle Eastern and Remaining groups are not included in this calculation (gnomAD v4.1 Joint)
* #### Allele frequency African/African American ((Allele_frequency_joint_AFR_GnomADv41))
Allele frequency in exomes and genomes of African/African American Populations (gnomAD v4.1 Joint)
* #### Allele frequency Admixed American ((Allele_frequency_joint_AMR_GnomADv41))
Allele frequency in exomes and genomes of Admixed American Populations (gnomAD v4.1 Joint)
* #### Allele frequency Amish ((Allele_frequency_joint_AMI_GnomADv41))
Allele frequency in genomes of Amish Populations. This group is only sequenced by genomes (gnomAD v4.1 Joint)
* #### Allele frequency Ashkenazi Jewish ((Allele_frequency_joint_ASJ_GnomADv41))
Allele frequency in exomes and genomes of Ashkenazi Jewish Populations (gnomAD v4.1 Joint)
* #### Allele frequency East Asian ((Allele_frequency_joint_EAS_GnomADv41))
Allele frequency in exomes and genomes of East Asian Populations (gnomAD v4.1 Joint)
* #### Allele frequency Finnish ((Allele_frequency_joint_FIN_GnomADv41))
Allele frequency in exomes and genomes of Finnish Populations (gnomAD v4.1 Joint)
* #### Allele frequency Middle Eastern ((Allele_frequency_joint_MID_GnomADv41))
Allele frequency in exomes and genomes of Middle Eastern Populations (gnomAD v4.1 Joint)
* #### Allele frequency Non-Finnish European ((Allele_frequency_joint_NFE_GnomADv41))
Allele frequency in exomes and genomes of Non-Finnish European Populations (gnomAD v4.1 Joint)
* #### Allele frequency South Asian ((Allele_frequency_joint_SAS_GnomADv41))
Allele frequency in exomes and genomes of South Asian Populations (gnomAD v4.1 Joint)
* #### Allele frequency Remaining ((Allele_frequency_joint_REMAINING_GnomADv41))
Allele frequency in exomes and genomes of individuals who were not assigned to any of the genetic ancestry groups above (gnomAD v4.1 Joint)

#### Graphical gnomAD Data
Graphical gnomAD data can be viewed by expanding any of the gnomAD (Graphical) subtiles. Two Graphs are available; one of the graphs is custom scaled to the allele frequencies by default (right side). Hovering over each bar will give you the numerical value represented in the population subset. You can click anywhere on the gnomAD (scaled) graph to change the scale between 1.0% (.01), 0.1% (.001), and the custom, default scale. This variety of scales will allow you to view all possible Allele Frequencies graphically.
Each group found on the x-axis of the bar chart can be found in the list of fields described in the gnomAD Populations section.

#### Numerical gnomAD Data
Numerical gnomAD data fields show numerical minor allele frequency data associated with each population, as well as the overall allele frequency. All of the minor allele frequencies are consistent with the graphs shown in the gnomAD (Graphical) subtile.

#### Exomic and Genomic Allele Frequencies for Variant Interpretation
There is not yet scientific consensus about whether exomic, genomic, or a combined allele frequency is most reliable for variant interpretation. The data display currently provides them separately. If you have expertise or opinions about how to use exomic and genomic allele frequencies together, please feel free to [contact us](mailto:brca-exchange-contact@genomicsandhealth.org?subject=BRCA%20Exchange%20Literature%20Search). 

#### gnomAD Flags
A flag on a variant can serve one of two purposes: either it indicates that there is something unusual about the data (such as low quality), or it indicates an annotation for the variant. For further information about why a variant might be flagged, use the gnomAD link in the tile to view the variant in gnomAD, along with detailed quality metrics. For more information about flags, visit the MacArthur Lab’s [blog post](https://macarthurlab.org/2018/10/17/gnomad-v2-1/) about the latest gnomAD release. 

#### Transcripts

The functional impact of a variant depends in part on where it appears in the RNA transcript.  The choice of transcript can affect the computation of the flags in gnomAD.  While gnomAD computes information relative to many transcripts, by default it shows information relative to the ENSEMBL canonical transcripts.  For _BRCA1_, this transcript is not clinically-relevant, and BRCA Exchange shows data from a different transcript.  For _BRCA2_, it shows data relative to gnomAD’s default transcript.  

The BRCA Exchange displays data associated with ENSEMBL transcript ENST00000357654 for BRCA1, and ENSMBL transcript ENST00000544455 for BRCA 2. These correspond to RefSeq transcript [NM_007294.3](http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_292.xml) and [NM_000059.3](http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_293.xml) respectively (per LRG).  Additional data, including detailed populations, quality scores, and flags relative to other transcripts, are available at gnomAD.

For more information about efforts to unify NCBI and EMBL-EBI datasets, please visit [MANE](https://ncbiinsights.ncbi.nlm.nih.gov/2019/03/12/mane-select-v0-5/) (Matched Annotation from the NCBI and EMBL-EBI). Though they have not yet covered _BRCA1_ and _BRCA2_, the MANE Select set will eventually comprise of “a matched representative transcript for every human protein-coding gene.”

#### GnomAD PopMax Filtering AF (95% Confidence) and Population
The PopMax Filtering Allele Frequency represents the 95% confidence threshold on the population frequencies, as evaluated over all of the populations.  This is usually the highest allele frequency for any population.  In cases where some population has a high observed allele frequency but low overall counts (i.e. low allele number), this observed allele frequency might not be statistically significant, and the filtering allele frequency might be the allele frequency of a different population.  The Popmax Filtering Allele Frequency population indicates the population from which the Popmax Filtering Allele Frequency was derived.  When studying a heritable disorder with a large disease cohort, if the observed allele frequency in the disease cohort exceeds the 95% filtering allele frequency, one can assume (with 95% confidence) that the variant in question is too common to be causative of the disease and can be dismissed.

##### References

1. [Whiffin et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/?term=28518168)

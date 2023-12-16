# Tremblay et al. 2024 custom code and data

This repository contains the custom code and data generated for the Tremblay et al. manuscript (in review) within the `preprocess/` and `analyses/` folders. The `data/` folder contains additional processed data not included in the supplementary of the published article. Please get in contact with the first author regarding any doubts or questions (`benjmtremblay [at] gmail [dot] com`). 

All bash scripts (i.e. scripts ending with `.sh`) can be run as: `bash script.sh`. All R scripts (i.e. scripts ending with `.R`) can be run as: `Rscript script.R`.

## `preprocess/`

This folder contains the scripts to transform the data from a set of raw FASTQ files to annotated peaks/transcripts/TSSs.

## `analyses/`

This folder contains R scripts for the generation of processed data analyzed in the manuscript.

## `data/`

This folder contains the main generated data in the course of the preprocessing and analysis steps.

## Computing environments, software used

All scripts in the `preprocess/` directory were run in a compute cluster with Intel(R) Xeon(R) Gold 6130 CPU, 1535 GB RAM, and running CentOS 7.7.1908 and require a few hours of runtime. All scripts in the `analysis/` directory were run using a MacBook Pro M1 with 16 GB RAM, and running macOS Big Sur 11.7.10 and require a few minutes of runtime.

The following software were used:

- [HOMER](http://homer.ucsd.edu/homer/index.html) (version 4.11.1)
- [STAR](https://github.com/alexdobin/STAR) (version 2.7.9a)
- [samtools](https://github.com/samtools/samtools) (version 1.14)
- [R](https://cran.r-project.org) (version 4.3.1)
- [GFF3sort](https://github.com/billzt/gff3sort) (version 1.0.0)
- [stringtie](https://ccb.jhu.edu/software/stringtie/) (version 2.1.4)
- [fastp](https://github.com/OpenGene/fastp) (version 0.23.4)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version 2.4.4)
- [MACS3](https://github.com/macs3-project/MACS) (version 3.0.0b3)
- [bedtools](https://github.com/arq5x/bedtools2) (version 2.30.0)
- [bigWigAverageOverBed](http://hgdownload.soe.ucsc.edu/admin/exe/) (version 2)
- [streme](https://meme-suite.org/meme/) (version 5.4.1)
- [tomtom](https://meme-suite.org/meme/) (version 5.4.1)
- [sea](https://meme-suite.org/meme/) (version 5.4.1)
- [Genrich](https://github.com/jsh58/Genrich) (version 0.6.1)
- [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) (version 0.12.6)
- [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml) (version 0.12.7)
- awk (version 20200816)
- python (version 3.11.4)

awk and python are generally preinstalled on most unix systems. Other
software must be installed one-by-one from their respective repositories.

The follow R packages were used:

- universalmotif (version 1.18.1)
- rtracklayer (version 1.60.0)
- Biostrings (version 2.69.2)
- csaw (version 1.34.0)
- edgeR (version 3.42.4)
- GenomicFeatures (version 1.52.1)
- readr (version 2.1.4)
- ggplot2 (version 3.4.2)
- matrixStats (version 1.0.0)
- WGCNA (version 1.72.1)
- BSgenome.Athaliana.TAIR.TAIR9 (version 1.3.1000)
- zoo (version 1.8.12)
- genomation (version 1.32.0)

To install these packages, install BiocManager with
`install.packages("BiocManager")`, then install packages with
`BiocManager::install("<PACKAGE_NAME>")`. This may require a few
minutes if compiled binaries are available, otherwise compiling
all packages from source may require over an hour of time.

## Additional data not included

Read density files and ACRs/TSSs/transcript files can be found in the GEO repository for this manuscript (*to be included upon publication*).

The following external data can be directly downloaded from the linked sources:

- [Araport11 GFF3 annotations](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release)
- [PhastCons scores](http://plantregmap.gao-lab.org/download.php#alignment-conservation)
- [PhyloP scores](http://plantregmap.gao-lab.org/download.php#alignment-conservation)
- [Leaf MNase-seq read density](https://bioinfor.yzu.edu.cn/download/plantdhs/Ath_leaf_NPS.bw)
- [Seedling Hi-C data](https://static-content.springer.com/esm/art%3A10.1038%2Fs41477-021-01004-x/MediaObjects/41477_2021_1004_MOESM6_ESM.xlsx)
- [Endosperm INT-Hi-C data](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/49/8/10.1093_nar_gkab191/1/gkab191_supplemental_files.zip?Expires=1688479597&Signature=x8HGjtgNXfv6alCNCqWDEzqNPAPr-jwTI7Ka7ncSay~J~Jf1nQ9947Jr0ikXsn4LAX0vsgsAS2ZOoFeF~DKAyI1VI3PPRIVqcQwZtXJkeTlER3IiDrnBJDb0ustA6SQN8IQL1~vjnCIbzVXhVOJwJj0vtrDi3xNxAaHpWhD2Hk-3yWsjrpBpOgZGsemJEoQCDXZker1SD0-Ubopu34neeHaAn2o07CpW2uko0MHmPCeE0cg9wtCVziWJpn0qG--TVCmY2DXO2EJ~LPOf~CQpYvU1TJan7TufamVs98eT-jbkaPqNh1FhgRPzBQN5mfgpM49sQHUTy1Mwf1NSLz8hNA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)
- [Seedling CAGE read density](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE136356&format=file)
- [ABI5 DAP-seq read density](http://systemsbiology.cau.edu.cn/chromstates/At_bwfile/ABI5-SRX670509.bw)

The following external datasets require additional processing such as read mapping and manual generation of bigWig files. Their sequence data repository accessions are as follows:

- Seedling total RNAPII ChIP-seq (Accession: DRA010413)
- Seedling RNAPII-Ser2P ChIP-seq (Accession: DRA010413)
- Seedling RNAPII-Ser5P ChIP-seq (Accession: DRA010413)
- Seedling H3K4me3 ChIP-seq (Accession: GSE96834)
- Seedling H3K9ac ChIP-seq (Accession: GSE79524)
- Seedling H3K4me1 ChIP-seq (Accession: DRA010413)
- Seedling H3K36me3 ChIP-seq (Accession: GSE96834)
- Seedling H2AZ ChIP-seq (Accession: GSE96834)
- Seedling H2AK121ub ChIP-seq (Accession: GSE89357)
- Seedling H3K27me3 ChIP-seq (Accession: GSE89357)
- Seedling GRO-cap (Accession: GSE83108)


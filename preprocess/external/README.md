# `preprocess/external/`

Processing of the external ChIP-seq, GRO-cap, GRO-seq and CAGE datasets used in the manuscript. The raw FASTQ files can be downloaded from SRA/DDBJ using the accessions listed in the script, and the CAGE bigWigs from GEO (GSE136356).

## Software used

- HOMER v4.11.1
- bowtie2 v2.4.4
- STAR v2.7.9a
- samtools v1.14
- BEDtools v2.30.0
- bigWigToBedGraph and bedGraphToBigWig (UCSC)

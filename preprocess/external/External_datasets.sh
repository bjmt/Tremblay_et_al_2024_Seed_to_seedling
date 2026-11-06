# Benjamin J.M. Tremblay
#
# Processing of external ChIP-seq, GRO-cap, GRO-seq and CAGE datasets

set -ex

mkdir -p external_tagDirs
mkdir -p external_bedGraphs
mkdir -p external_bigWigs

# ChIP-seq ---------------------------------------------------------------------
#
# RNAPII (DRA010413): DRR235321, DRR235325
# RNAPII-Ser2P (DRA010413): DRR235329, DRR235331
# RNAPII-Ser5P (DRA010413): DRR235330, DRR235332
# H3K4me1 (DRA010413): DRR235274
# H3K4me3 (GSE96834): SRR5364425
# H3K36me3 (GSE96834): SRR5364427, SRR5364428
# H2AZ (GSE96834): SRR5364421, SRR5364422
# H3K9ac (GSE79524): SRR5224475
# H3K27me3 (GSE89357): SRR4734660, SRR4734662, SRR5278091, SRR5278093
# H2AK121ub (GSE89357): SRR4734656, SRR4734658, SRR4734664, SRR4734666,
#   SRR4734668, SRR4734670

CHIP_RUNS=("DRR235321" "DRR235325" "DRR235329" "DRR235330" "DRR235331"
  "DRR235332" "DRR235274" "SRR5364425" "SRR5364427" "SRR5364428" "SRR5364421"
  "SRR5364422" "SRR5224475" "SRR4734660" "SRR4734662" "SRR5278091" "SRR5278093"
  "SRR4734656" "SRR4734658" "SRR4734664" "SRR4734666" "SRR4734668"
  "SRR4734670")

for i in ${CHIP_RUNS[@]} ; do

  # Trim reads with HOMER

  gzip -cd ${i}.fastq.gz > ${i}.fastq

  homerTools trim \
    -3 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -mis 2 \
    -minMatchLength 4 \
    -min 20 \
    -suffix trimmed.fq \
    ./${i}.fastq

  rm ${i}.fastq

  # Map with bowtie2

  bowtie2 -x tair10 ./${i}.fastq.trimmed.fq > ${i}.sam

  gzip ${i}.fastq.trimmed.fq

  # Generate tagDirs with HOMER

  makeTagDirectory external_tagDirs/${i}-tagDir ${i}.sam -genome tair10

  samtools view -b ${i}.sam > ${i}.bam
  rm ${i}.sam

done

# Merge replicates

makeTagDirectory external_tagDirs/RNAPII_14D_merged-tagDir \
  -d external_tagDirs/DRR235321-tagDir external_tagDirs/DRR235325-tagDir
makeTagDirectory external_tagDirs/RNAPIISer2P_14D_merged-tagDir \
  -d external_tagDirs/DRR235329-tagDir external_tagDirs/DRR235331-tagDir
makeTagDirectory external_tagDirs/RNAPIISer5P_14D_merged-tagDir \
  -d external_tagDirs/DRR235330-tagDir external_tagDirs/DRR235332-tagDir
makeTagDirectory external_tagDirs/H3K36me3_10D_merged-tagDir \
  -d external_tagDirs/SRR5364427-tagDir external_tagDirs/SRR5364428-tagDir
makeTagDirectory external_tagDirs/H2AZ_10D_merged-tagDir \
  -d external_tagDirs/SRR5364421-tagDir external_tagDirs/SRR5364422-tagDir
makeTagDirectory external_tagDirs/H3K27me3_7D_merged-tagDir \
  -d external_tagDirs/SRR4734660-tagDir external_tagDirs/SRR4734662-tagDir \
     external_tagDirs/SRR5278091-tagDir external_tagDirs/SRR5278093-tagDir
makeTagDirectory external_tagDirs/H2AK121ub_7D_merged-tagDir \
  -d external_tagDirs/SRR4734656-tagDir external_tagDirs/SRR4734658-tagDir \
     external_tagDirs/SRR4734664-tagDir external_tagDirs/SRR4734666-tagDir \
     external_tagDirs/SRR4734668-tagDir external_tagDirs/SRR4734670-tagDir

mv external_tagDirs/DRR235274-tagDir external_tagDirs/H3K4me1_14D-tagDir
mv external_tagDirs/SRR5364425-tagDir external_tagDirs/H3K4me3_10D-tagDir
mv external_tagDirs/SRR5224475-tagDir external_tagDirs/H3K9ac_12D-tagDir

# GRO-cap and GRO-seq (GSE83108) ------------------------------------------------
#
# GRO-cap: SRR3647033
# GRO-seq: SRR3647034, SRR3647035

GRO_RUNS=("SRR3647033" "SRR3647034" "SRR3647035")

for i in ${GRO_RUNS[@]} ; do

  homerTools trim \
    -3 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -mis 2 \
    -minMatchLength 4 \
    -min 20 \
    ${i}.fastq.gz

  STAR --genomeDir STAR-tair10 --runThreadN 6 \
    --readFilesIn ${i}.fastq.gz.trimmed \
    --outFileNamePrefix ${i}. \
    --outSAMstrandField intronMotif \
    --outMultimapperOrder Random \
    --outSAMmultNmax 1 \
    --outFilterMultimapNmax 10000 \
    --limitOutSAMoneReadBytes 10000000

  samtools sort -o ${i}.sorted.bam ${i}.Aligned.out.sam

done

makeTagDirectory external_tagDirs/GROcap-tagDir SRR3647033.Aligned.out.sam \
  -genome tair10 -checkGC -fragLength 150
makeTagDirectory external_tagDirs/GROseq1-tagDir SRR3647034.Aligned.out.sam \
  -genome tair10 -checkGC -fragLength 150
makeTagDirectory external_tagDirs/GROseq2-tagDir SRR3647035.Aligned.out.sam \
  -genome tair10 -checkGC -fragLength 150

samtools merge GROseq.merged.sorted.bam SRR3647034.sorted.bam SRR3647035.sorted.bam

STAR --runMode inputAlignmentsFromBAM \
  --inputBAMfile SRR3647033.sorted.bam \
  --outWigType bedGraph read1_5p \
  --outWigStrand Stranded \
  --outWigNorm RPM \
  --outFileNamePrefix external_bedGraphs/GROcap.5p.RPM.

STAR --runMode inputAlignmentsFromBAM \
  --inputBAMfile GROseq.merged.sorted.bam \
  --outWigType bedGraph \
  --outWigStrand Stranded \
  --outWigNorm RPM \
  --outFileNamePrefix external_bedGraphs/GROseq.merged.RPM.

for j in external_bedGraphs/*.bg ; do
  mv $j ${j%.bg}.bedGraph
done

gzip external_bedGraphs/*.bedGraph

# CAGE (GSE136356) ---------------------------------------------------------------

CAGE_SAMPLES=("wt" "hen2" "rrp4")
CAGE_REPS=(1 2 3)

for i in ${CAGE_SAMPLES[@]} ; do

  for j in ${CAGE_REPS[@]} ; do

    bigWigToBedGraph *${i}_0_R${j}.raw.plus.bw tmp_plus.bedGraph
    bigWigToBedGraph *${i}_0_R${j}.raw.minus.bw tmp_minus.bedGraph

    awk 'BEGIN{IFS=OFS="\t"}{print $1,$2,$3,".",$4,"+"}' \
      tmp_plus.bedGraph > tmp.bed
    awk 'BEGIN{IFS=OFS="\t"}{print $1,$2,$3,".",$4,"-"}' \
      tmp_minus.bedGraph >> tmp.bed

    bedtools sort -i tmp.bed > tmp_sorted.bed

    rm tmp_plus.bedGraph tmp_minus.bedGraph tmp.bed

    makeTagDirectory external_tagDirs/${i}_${j}.cage-tagDir \
      tmp_sorted.bed -format bed -force5th

    rm tmp_sorted.bed

  done

done

makeTagDirectory external_tagDirs/CAGE_14D_merged-tagDir \
  -d external_tagDirs/wt_1.cage-tagDir external_tagDirs/wt_2.cage-tagDir \
     external_tagDirs/wt_3.cage-tagDir

# bedGraphs and bigWigs --------------------------------------------------------

FINAL_TAGDIRS=("RNAPII_14D_merged" "RNAPIISer2P_14D_merged"
  "RNAPIISer5P_14D_merged" "H3K4me1_14D" "H3K4me3_10D" "H3K9ac_12D"
  "H3K36me3_10D_merged" "H2AZ_10D_merged" "H3K27me3_7D_merged"
  "H2AK121ub_7D_merged" "CAGE_14D_merged")

for i in ${FINAL_TAGDIRS[@]} ; do

  makeUCSCfile external_tagDirs/${i}-tagDir > external_bedGraphs/${i}.bedGraph
  gzip external_bedGraphs/${i}.bedGraph

  gzip -cd external_bedGraphs/${i}.bedGraph.gz \
    | bedtools sort -g chrom.sizes -i stdin \
    > tmp.bedGraph
  bedGraphToBigWig tmp.bedGraph chrom.sizes external_bigWigs/${i%_merged}.bw
  rm tmp.bedGraph

done

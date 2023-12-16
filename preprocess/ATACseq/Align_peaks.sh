# Benjamin J.M. Tremblay, June 2023
#
# Trimming, mapping, filtering, peak calling and quantification of ATAC-seq data

SAMPLES=("DS1" "DS2" "L61" "L62" "L261" "L262" "L571" "L572")
GROUPS=("DS" "L6" "L26" "L57")

mkdir -p atacseq_trimmed
mkdir -p atacseq_bams
mkdir -p atacseq_genrich
mkdir -p atacseq_bed
mkdir -p atacseq_peaks
mkdir -p atacseq_bed1
mkdir -p atacseq_bedGraphs
mkdir -p atacseq_counts
mkdir -p atacseq_bw

for i in ${SAMPLES[@]} ; do

  # Trim reads with fastp

  fastp \
    -a CTGTCTCTTATACACATCT \
    -i atacseq_fastq/${i}*R1_001.fastq.gz \
    -I atacseq_fastq/${i}*R2_001.fastq.gz \
    -o atacseq_trimmed/${i}_1.trimmed.fastq.gz \
    -O atacseq_trimmed/${i}_2.trimmed.fastq.gz \
    -j atacseq_trimmed/${i}_fastp.json \
    -h atacseq_trimmed/${i}_fastp.html

  # Map reads with bowtie2, sort with samtools

  bowtie2 --very-sensitive -X 2000 --dovetail \
    -x bowtie2-tair10/tair10 \
    -1 atacseq_trimmed/${i}_1.trimmed.fastq.gz \
    -2 atacseq_trimmed/${i}_2.trimmed.fastq.gz \
    | samtools view -u - \
    | samtools sort -n -o atacseq_bams/${i}.bam

  # Filter nuclear reads only and MAPQ of at least 10 with samtools

  samtools view -q 10 -L nuclear.bed -hbo atacseq_bams/${i}.nuclear.bam \
    atacseq_bams/${i}.bam

done

for i in ${SAMPLES[@]} ; do

  # Use Genrich to identify duplicates

  Genrich -t rnaseq_bams/${i}.nuclear.bam \
    -o atacseq_genrich/${i}.narrowPeak \
    -f atacseq_genrich/${i}.log \
    -b atacseq_genrich/${i}.frags.bed \
    -k atacseq_genrich/${i}.pileup.bed \
    -R atacseq_genrich/${i}.duplicates.txt \
    -r -e Pt,Mt -E blacklist.bed \
    -m 10 -j -q 0.05 -z -v -a 0 -y

done

for i in ${SAMPLES[@]} ; do

  # Use samtools and getReads.py from Genrich to remove duplicates from BAM

  samtools view -h atacseq_bams/${i}.nuclear.bam \
    | python getReads.py - no atacseq_genrich/${i}.duplicates.txt.gz - \
    | samtools view -b - > atacseq_bams/${i}.nuclear.dedup.bam

  # Sort with samtools

  samtools sort -o atacseq_bams/${i}.nuclear.dedup.sort.bam \
    atacseq_bams/${i}.nuclear.dedup.bam

  rm atacseq_bams/${i}.nuclear.dedup.bam

  # Use BEDtools to filter reads from blacklist regions and generate a BED file

  bedtools bamtobed -i atacseq_bams/${i}.nuclear.dedup.sort.bam \
    | bedtools subtract -A -a stdin -b blacklist.bed \
    | gzip > atacseq_bed/${i}.nuclear.dedup.sort.bed.gz

done

for i in ${SAMPLES[@]} ; do

  # Call peaks with MACS

  macs3 callpeak --nomodel --shift -100 --extsize 200 \
    -t atacseq_bed/${i}.nuclear.dedup.sort.bed.gz \
    -n ${i} -g 1.191e8 --keep-dup all \
    --outdir atacseq_peaks -p 0.05 -f BED --call-summits

done

for i in ${GROUPS[@]} ; do

  # Merge peaks between replicates and keep those present in both with BEDtools

  PEAK1=atacseq_peaks/${i}1_peaks.narrowPeak
  PEAK1=atacseq_peaks/${i}2_peaks.narrowPeak

  awk '{if($9>=1.30103) print $0}' $PEAK1 > tmp_merged.narrowPeak
  awk '{if($9>=1.30103) print $0}' $PEAK2 >> tmp_merged.narrowPeak

  bedtools sort -i tmp_merged.narrowPeak \
    | bedtools merge -i stdin -c 4 -o collapse \
    | grep ${i}1 | grep ${i}2 \
    > atacseq_peaks/${i}.merged.narrowPeak

done

rm tmp_merged.narrowPeak

# Get final peak set from all samples with BEDtools

cat atacseq_peaks/*.merged.narrowPeak \
  | cut -f1-3 \
  | bedtools sort -i stdin \
  | bedtools merge -i stdin \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"ATAC"$1"_"int(($2+$3)/2)}' \
  > atacseq_peaks/peaks_unadjusted.bed

# Get read 5' ends with BEDtools

for i in ${SAMPLES[@]} ; do

  bedtools bamtobed -i atacseq_bams/${i}.nuclear.dedup.sort.bam \
    | awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else if($6=="-"){print $1,$3-1,$3,$4,$5,$6}}' \
    | gzip > atacseq_bed1/${i}.size1.bed.gz

done

# Merge all samples with BEDtools

gzcat atacseq_bed1/*.size1.bed.gz \
  | bedtools sort -i stdin \
  | gzip > atacseq_bed1/all_merged.size1.bed.gz

# BED->bedGraph with BEDtools

gzcat atacseq_bed1/all_merged.size1.bed.gz \
  | bedtools slop -l 100 -r 99 -s -g chrom.sizes -i stdin \
  | bedtools genomecov -i stdin -bg -g chrom.sizes \
  | gzip > atacseq_bedGraphs/all_merged.bedGraph.gz

# Get peak summits, rename peaks and make a version of peaks 500 bp wide in R

Rscript Align_peaks_helper_part1.R

# Make a thick BED

awk 'BEGIN{IFS=OFS="\t"}{s1=$4;int(sub(/ATAC[12345]_/,"",s1));print $1,$2,$3,$4,1,".",s1-1,S1}' \
  atacseq_peaks/peaks_final.bed > atacseq_peaks/peaks_final_thick.bed

# Get peak counts per sample

for i in ${SAMPLES[@]} ; do

  zcat atacseq_bed1/${i}.size1.bed.gz \
    | bedtools coverage -counts -a atacseq_peaks/peaks_final.bed -b stdin \
    | cut -f4-5 \
    > atacseq_counts/${i}.counts.txt

done

# Merge peak counts

paste \
  <(awk 'BEGIN{print "PeakID\tDS1"}{print}' atacseq_counts/DS1.counts.txt) \
  <(cut -f2 atacseq_counts/DS2.counts.txt | awk 'BEGIN{print "DS2"}{print}') \
  <(cut -f2 atacseq_counts/L61.counts.txt | awk 'BEGIN{print "L61"}{print}') \
  <(cut -f2 atacseq_counts/L62.counts.txt | awk 'BEGIN{print "L62"}{print}') \
  <(cut -f2 atacseq_counts/L261.counts.txt | awk 'BEGIN{print "L261"}{print}') \
  <(cut -f2 atacseq_counts/L262.counts.txt | awk 'BEGIN{print "L262"}{print}') \
  <(cut -f2 atacseq_counts/L571.counts.txt | awk 'BEGIN{print "L571"}{print}') \
  <(cut -f2 atacseq_counts/L572.counts.txt | awk 'BEGIN{print "L572"}{print}') \
  > atacseq_counts/all.counts.txt

# Tile the genome in 10 bp bins

bedtools makewindows -g chrom.sizes -w 10 > genome_tiled10bp.bed

# Merge replicate coverage in 10 bp bins

for i in ${GROUPS[@]} ; do

  zcat atacseq_bed1/${i}[12].size1.bed.gz \
    | bedtools sort -i stdin \
    | bedtools slop -i stdin -g chrom.sizes -l 25 -r 24 \
    | bedtools coverage -counts -a genome_tiled10bp.bed -b stdin \
    | gzip > atacseq_bedGraphs/${i}.tiled10bp.bedGraph.gz

done

# Get normalized bigWigs

Rscript Align_peaks_helper_part2.R


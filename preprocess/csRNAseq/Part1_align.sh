# Benjamin J.M. Tremblay, June 2023
#
# Trimming, mapping, Tag directory and bedGraph creation for the (cs)RNA-seq

SAMPLES=("DS1" "DS2" "S241" "S242" "S721" "S722" "L61" "L62" "L261" "L262" "L571" "L572" "hen21" "hen22" "rrp41" "rrp42")
GROUPS=("DS" "S24" "S72" "L6" "L26" "L57" "hen2" "rrp4")

mkdir -p srnaseq_trimmed
mkdir -p srnaseq_bams
mkdir -p srnaseq_bams_filtered
mkdir -p srnaseq_bedGraphs

mkdir -p csrnaseq_trimmed
mkdir -p csrnaseq_bams
mkdir -p csrnaseq_bams_filtered
mkdir -p csrnaseq_bedGraphs

mkdir -p rnaseq_trimmed
mkdir -p rnaseq_bams
mkdir -p rnaseq_bedGraphs
mkdir -p rnaseq_tagDirs

# csRNA-seq / sRNA-seq:

for i in ${SAMPLES[@]} ; do

  for j in srnaseq csrnaseq ; do
 
    # Trim raw reads with HOMER

    homerTools trim \
      -3 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      -mis 2 \
      -minMatchLength 4 \
      -min 20 \
      ${j}_fastq/${i}*R1_001.fastq.gz

    gzip ${j}_fastq/${i}*trimmed

    mv ${j}_fastq/${i}*R1_001.fastq.gz.trimmed.gz \
      ${j}_trimmed/${i}_R1_001.trimmed.fastq.gz

  done

done

for i in ${SAMPLES[@]} ; do

  for j in srnaseq csrnaseq ; do
 
    # Map to the genome with STAR

    STAR \
      --genomeDir STAR-tair10 \
      --readFilesCommand zcat \
      --readFilesIn ${j}_trimmed/${i}_R1_001.trimmed.fastq.gz \
      --outFileNamePrefix ${j}_bams/${i}. \
      --outSAMstrandField intronMotif \
      --outMultimapperOrder Random \
      --outSAMmultNmax 1 \
      --outFilterMultimapNmax 10000 \
      --limitOutSAMoneReadBytes 10000000

    # Sort and convert sam->bam with samtools

    samtools view -bhS ${j}_bams/${i}*.sam \
      | samtools sort - \
      > ${j}_bams/${i}.sorted.bam

    rm ${j}_bams/${i}*.sam

    # Filter reads
    # - Max 1 mismatch
    # - MAPQ at least 10
    # - No gaps
    # - No soft clipping
    # - Read length from 20 to 70 only

    samtools view -t chrom.sizes ${j}_bams/${i}.sorted.bam \
      | grep "nM:i:[01]" \
      | awk '$5>=10{print $0}' \
      | awk '$6!~/N/{print $0}' \
      | awk '$6!~/S/{print $0}' \
      | awk 'length($10)>=20||length($10)<=70{print $0}' \
      > ${j}_bams_filtered/${i}.filtered.sam

    # Generate tagDirs with HOMER

    makeTagDirectory ${j}_tagDirs/${i}-tagDir/ \
      ${j}_bams_filtered/${i}.filtered.sam \
      -genome tair10 \
      -checkGC

    # Convert sam->bam with samtools

    samtools view -b -t chrom.sizes ${j}_bams_filtered/${i}.filtered.sam \
      > ${j}_bams_filtered/${i}.filtered.bam

    rm ${j}_bams_filtered/${i}.filtered.sam

    # Generate bedGraphs with HOMER

    makeUCSCfile ${j}_tagDirs/${i}-tagDir \
      -style tss -strand + > ${j}_bedGraphs/${i}.plus.bedGraph
    gzip -f ${j}_bedGraphs/${i}.plus.bedGraph

    makeUCSCfile ${j}_tagDirs/${i}-tagDir \
      -style tss -strand - > ${j}_bedGraphs/${i}.minus.bedGraph
    gzip -f ${j}_bedGraphs/${i}.minus.bedGraph

  done

done

# Get replicate-merged bedGraphs:

for i in ${GROUPS[@]} ; do

  samtools merge csrnaseq_bams_filtered/${i}.merged.bam \
    csrnaseq_bams_filtered/${i}[12].filtered.bam

  makeTagDirectory csrnaseq_tagDirs/${i}.merged-tagDir/ \
    csrnaseq_bams_filtered/${i}.merged.bam \
    -genome tair10

  makeUCSCfile csrnaseq_tagDirs/${i}.merged-tagDir/ \
    -style tss -strand + > csrnaseq_bedGraphs/${i}.merged.plus.bedGraph
  gzip -f csrnaseq_bedGraphs/${i}.merged.plus.bedGraph

  makeUCSCfile csrnaseq_tagDirs/${i}.merged-tagDir/ \
    -style tss -strand - > csrnaseq_bedGraphs/${i}.merged.minus.bedGraph
  gzip -f csrnaseq_bedGraphs/${i}.merged.minus.bedGraph

done

# RNA-seq:

for i in ${SAMPLES[@]} ; do

  # Trim raw reads with HOMER

  homerTools trim
    -3 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -mis 2 \
    -minMatchLength 4 \
    -min 20 \
    -pe \
    rnaseq_fastq/${i}*R1_001.fastq.gz \
    rnaseq_fastq/${i}*R2_001.fastq.gz

  gzip rnaseq_fastq/${i}*trimmed

  mv rnaseq_fastq/${i}*R1_001.fastq.gz.trimmed.gz \
    rnaseq_trimmed/${i}_R1_001.trimmed.fastq.gz
  mv rnaseq_fastq/${i}*R2_001.fastq.gz.trimmed.gz \
    rnaseq_trimmed/${i}_R2_001.trimmed.fastq.gz

done

for i in ${SAMPLES[@]} ; do

  # Map to the genome with STAR

  STAR \
    --genomeDir STAR-tair10 \
    --readFilesCommand zcat \
    --readFilesIn rnaseq_trimmed/${i}_R[12]_001.trimmed.fastq.gz \
    --outFileNamePrefix rnaseq_bams/${i}. \
    --outSAMstrandField intronMotif \
    --outMultimapperOrder Random \
    --outSAMmultNmax 1 \
    --outFilterMultimapNmax 10000 \
    --limitOutSAMoneReadBytes 10000000

  # Generate tagDirs with HOMER

  makeTagDirectory rnaseq_tagDirs/${i}_r2-tagDir/ \
    rnaseq_bams/${i}*.sam \
    -genome tair10 \
    -read2 \
    -checkGC

  # Sort and convert sam->bam with samtools

  samtools view -bhS rnaseq_bams/${i}*.sam \
    | samtools sort - \
    > rnaseq_bams/${i}.sorted.bam

  rm rnaseq_bams/${i}*.sam

  # Generate bedGraphs with STAR

  STAR \
    --runMode inputAlignmentsFromBAM \
    --inputBAMfile rnaseq_bams/${i}.sorted.bam \
    --outWigType bedGraph \
    --outWigNorm RPM \
    --outFileNamePrefix rnaseq_bedGraphs/${i}.

done

# Get replicate-merged bedGraphs:

for i in ${GROUPS[@]} ; do

  samtools merge rnaseq_bams/${i}.merged.bam \
    rnaseq_bams/${i}[12].sorted.bam

  STAR \
    --runMode inputAlignmentsFromBAM \
    --inputBAMfile rnaseq_bams/${i}.merged.bam \
    --outWigType bedGraph \
    --outWigNorm RPM \
    --outFileNamePrefix rnaseq_bedGraphs/${i}.merged.

done


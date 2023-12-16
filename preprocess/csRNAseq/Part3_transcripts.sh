# Benjamin J.M. Tremblay, June 2023
#
# De novo transcript identification and quantification

GROUPS=("DS" "S24" "S72" "L6" "L26" "L57" "hen2" "rrp4")

mkdir rnaseq_denovo
mkdir rnaseq_quantification

for i in ${GROUPS[@]} ; do

  # Merge replicate RNA-seq BAMs with samtools

  samtools merge rnaseq_bams/${i}.merged.bam \
    rnaseq_bams/${i}[12].bam

  # Filter merged BAM with samtools
  # - MAPQ at least 10
  # - No gaps
  # - No soft clipping

  samtools view -t chrom.sizes rnaseq_bams/${i}.merged.bam \
    | awk '$5>=10{print $0}' \
    | awk '$6!~/N/{print $0}' \
    | awk '$6!~/S/{print $0}' \
    > rnaseq_bams/${i}.merged.filtered.sam

  # Make tag directory with HOMER

  makeTagDirectory rnaseq_tagDirs/${i}_denovo-tagDir \
    rnaseq_bams/${i}.merged.filtered.sam -sspe

  # Convert sam->bam with samtools

  samtools view -t chrom.sizes -b rnaseq_bams/${i}.merged.filtered.sam \
    > rnaseq_bams/${i}.merged.filtered.bam

  rm rnaseq_bams/${i}.merged.filtered.sam

  # Find new transcripts with HOMER

  findPeaks rnaseq_tagDirs/${i}_denovo-tagDir \
    -style groseq -tssSize 100 -minBodySize 100 -endFold 100 -rev \
    -confPvalue 0.001 \
    -gtf rnaseq_denovo/${i}.homer.gtf

done

# Merge HOMER transcripts with gffcompare, remove overlaps with existing annotations

gffcompare -r Araport11_GFF3_genes_transposons.gff3 \
  rnaseq_denovo/*.homer.gtf \
  -o rnaseq_denovo/homer.merged -p HOMER_RNA

# Get primary TSS coordinates ready for stringtie

awk '{print $1,$3,$6,"TSS"}' csrnaseq_TSSs_final/tss.size1.bed \
  > rnaseq_denovo/tss.pos

for i in ${GROUPS[@]} ;

  # Find new transcripts with StringTie

  stringtie rnaseq_bams/${i}.merged.bam \
    -o rnaseq_denovo/${i}.stringtie.gtf \
    -G Araport11_GFF3_genes_transposons.gff3 \
    --rf -m 150 -s 1 -l ${i}.STRG -f 0.01 \
    --ptf rnaseq_denovo/tss.pos \
    -C rnaseq_denovo/${i}.stringtie_cov_refs.gtf \
    -A rnaseq_denovo/${i}.stringtie_gene_abund.txt

do

# merge StringTie transcripts with gffcompare, remove overlaps with existing annotations

gffcompare -r Araport11_GFF3_genes_transposons.gff3 \
  rnaseq_denovo/*.stringtie.gtf \
  -o rnaseq_denovo/stringtie.merged -p STRINGTIE_RNA

# Combine HOMER and StringTie transcripts in R

Rscript Part3_transcripts_helper1.R
Rscript Part3_transcripts_helper2.R

# GTF->GFF with gffread

gffread -E -O rnaseq_denovo/novel.transcripts.filtered.gtf \
  -o rnaseq_denovo/novel.transcripts.filtered.gff3

# GFF3->BED

< rnaseq_denovo/novel.transcripts.filtered.gff3 \
  awk '$0!~/^#/' \
  | awk 'BEGIN{OFS="\t"}/\ttranscripts\t/{print $1,$4-1,$5,$9,".",$7}' \
  | sed 's/ID=//g' \
  > rnaseq_denovo/novel.transcripts.filtered.bed

# Merge with Araport11 annotations and sort with GFF3sort

cp -f Araport11_GFF3_genes_transposons.gff3 Araport11_novel.unsorted.gff3
grep "^[^#]" rnaseq_denovo/novel.transcripts.filtered.gff3 \
  >> Araport11_novel.unsorted.gff3
gff3sort.pl --precise Araport11_novel.unsorted.gff3 \
  > Araport11_novel.gff3
rm Araport11_novel.unsorted.gff3

# Quantify all transcripts with StringTie

for i in ${GROUPS[@]} ; do
  for j in 1 2 ; do

    stringtie rnaseq_bams/${i}${j}.sorted.bam \
      -o rnaseq_quantification/${i}${j}.gtf \
      -G Araport11_novel.gff3 --rf -e \
      -A rnaseq_quantification/${i}${j}.gene_abund.txt

    printf "${i}${j} rnaseq_quantification/${i}${j}.gtf\n" \
      >> rnaseq_quantification/gtfs.txt

  done
done

# Merge all quantifications with prepDE.py from StringTie

prepDE.py -v -l 125 -i rnaseq_quantification/gtfs.txt \
  -g rnaseq_quantification/gene_counts.csv \
  -t rnaseq_quantification/transcript_counts.csv

# Convert counts csv->tsv

tr ',' '\t' < rnaseq_quantification/gene_counts.csv \
  > rnaseq_quantification/gene_counts.tsv
tr ',' '\t' < rnaseq_quantification/transcript_counts.csv \
  > rnaseq_quantification/transcript_counts.tsv


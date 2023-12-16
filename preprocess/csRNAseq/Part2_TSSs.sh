# Benjamin J.M. Tremblay, June 2023
#
# Calling TSS peaks

SAMPLES=("DS1" "DS2" "S241" "S242" "S721" "S722" "L61" "L62" "L261" "L262" "L571" "L572" "hen21" "hen22" "rrp41" "rrp42")
GROUPS=("DS" "S24" "S72" "L6" "L26" "L57" "hen2" "rrp4")

mkdir -p csrnaseq_TSSs
mkdir -p csrnaseq_TSSs_final

for i in ${SAMPLES[@]} ; do

  # Find TSSs with HOMER

  findcsRNATSS.pl \
    csrnaseq_tagDirs/${i}-tagDir \
    -includeSingleExons \
    -genome tair10 \
    -ntagThreshold 1 -noFilterRNA -minDistDiff 0.01 -L 1.5 \
    -i srnaseq_tagDirs/${i}-tagDir \
    -rna rnaseq_tagDirs/${i}_r2-tagDir \
    -o csrnaseq_TSSs/${i}

  # Convert to BED

  awk 'BEGIN{IFS="\t";OFS="\t"}NR>1{print $2,$3,$4,$1,$6,$5}' \
    csrnaseq_TSSs/${i}.tss.txt > csrnaseq_TSSs/${i}.tss.bed

done

for i in ${GROUPS[@]} ; do

  # Merge peaks between replicates with HOMER, only keep those present in both

  mergePeaks -strand csrnaseq_TSSs/${i}[12].tss.txt \
    | awk 'NR==1{print $0}(NR>1)&&($8==2){print $0}' \
    > csrnaseq_TSSs/${i}.tss.merged.txt

  # Get merged bedGraphs with bedtools

  zcat csrnaseq_bedGraphs/${i}[12].plus.bedGraph.gz \
    | sortBed -i stdin \
    | mergeBed -c 4 -o sum -i stdin \
    | gzip \
    > csrnaseq_bedGraphs/${i}.merged.plus.bedGraph.gz

  zcat csrnaseq_bedGraphs/${i}[12].minus.bedGraph.gz \
    | sortBed -i stdin \
    | mergeBed -c 4 -o sum -i stdin \
    | gzip \
    > csrnaseq_bedGraphs/${i}.merged.minus.bedGraph.gz

done

# Assemble final list of TSSs from all samples with HOMER, only keep nuclear

mergePeaks -strand \
  csrnaseq_TSSs/*.tss.merged.txt \
  | awk 'NR==1{print $0}NR>1{if($2!="Mt"&&$2!="Pt"){print $0}}' \
  > csrnaseq_TSSs_final/all.tss.merged.filtered.txt

# Convert to BED format with HOMER, sort with bedtools

pos2bed.pl -float \
  csrnaseq_TSSs_final/all.tss.merged.filtered.txt \
  | bedTools sort -i stdin \
  > csrnaseq_TSSs_final/all.tss.merged.filtered.bed

# Get raw read counts with HOMER

annotatePeaks.pl \
  csrnaseq_TSSs_final/all.tss.merged.filtered.txt \
  tair10 -gff3 Araport11_GFF3_genes_transposons.gff3 \
  -strand + -fragLength 1 -raw -d csrnaseq_tagDirs/*tagDir \
  > csrnaseq_TSSs_final/all.tss.merged.counts.txt

# Merged all bedGraphs into one with bedtools

zcat csrnaseq_bedGraphs/*merged.plus.bedGraph.gz \
  | sortBed -i stdin \
  | mergeBed -c 4 -o sum -i stdin \
  | gzip \
  > csrnaseq_bedGraphs/all.plus.bedGraph.gz

zcat csrnaseq_bedGraphs/*merged.minus.bedGraph.gz \
  | sortBed -i stdin \
  | mergeBed -c 4 -o sum -i stdin \
  | gzip \
  > csrnaseq_bedGraphs/all.minus.bedGraph.gz

# Use the bedGraphs merged from all samples to refine TSSs in R

Rscript Part2_helper.R

# Use the refined TSSs and the summits to create a thick BED

awk 'BEGIN{IFS=OFS="\t"}{s1=$4;int(sub(/TSS[12345]_/,"",s1));print $0,s1-1,s1}' \
  csrnaseq_TSSs_final/tss.shrunk.bed > csrnaseq_TSSs_final/tss.shrunk_size1.bed


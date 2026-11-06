# Benjamin J.M. Tremblay
#
# Motif tables for the figures (MEME suite 5.4.1, and minimotif and its
# to_bed.sh script from https://github.com/bjmt/yamtk for the genome scan):
# bash Part05_motifs_helper3.sh streme|sea|scan

set -exo pipefail

if [ "$1" == "streme" ] ; then

  mkdir -p streme_top2k

  for i in A1 A2 A3 A4 A5 ; do
    streme \
      --oc streme_top2k/${i} \
      --maxw 12 \
      --seed 1234 \
      --align center \
      --n fasta2/atac_bkg.fa \
      --p fasta2/atac_${i}.fa
    centrimo \
      --oc streme_top2k/${i}_centrimo \
      --neg fasta2/atac_bkg.fa \
      --flip \
      --optimize-score \
      fasta2/atac_${i}.fa \
      streme_top2k/${i}/streme.txt
  done

  for i in C1 C2 C3 C4 C5 C6 ; do
    streme \
      --oc streme_top2k/${i} \
      --maxw 12 \
      --seed 1234 \
      --align right \
      --n fasta2/tss_bkg.fa \
      --p fasta2/tss_${i}.fa
    centrimo \
      --oc streme_top2k/${i}_centrimo \
      --neg fasta2/tss_bkg.fa \
      --local \
      --optimize-score \
      fasta2/tss_${i}.fa \
      streme_top2k/${i}/streme.txt
  done

elif [ "$1" == "sea" ] ; then

  mkdir -p streme_final

  for i in A1 A2 A3 A4 A5 ; do
    sea \
      --oc streme_final/all_${i}_sea \
      --align center \
      --m all_final.meme \
      --p fasta2/atac_${i}.fa \
      --n fasta2/atac_bkg.fa
  done

  for i in C1 C2 C3 C4 C5 C6 ; do
    sea \
      --oc streme_final/all_${i}_sea \
      --align center \
      --m all_final.meme \
      --p fasta2/tss_${i}.fa \
      --n fasta2/tss_bkg.fa
  done

  tomtom \
    -oc streme_final/all_final_tomtom \
    -internal \
    all_final.meme \
    PlantTFDB_Ath_TF_binding_motifs.meme

elif [ "$1" == "scan" ] ; then

  minimotif \
    -m all_final.meme \
    -s genome.fa \
    -t 0.0001 -g -v \
  | to_bed.sh \
  | gzip > final_streme_motifs.bed.gz

fi


# Benjamin J.M. Tremblay, June 2023
#
# Motif analyses

set -ex

for i in C1 C2 C3 C4 C5 C6 ; do

  sea \
    --oc sea/sea_${i} \
    --align center \
    --m Final_motifs.meme \
    --p Promoters/Promoters_${i}.fa \
    --n Promoters/Promoters_bkg.fa

done

for i in A1 A2 A3 A4 A5 ; do

  sea \
    --oc sea/sea_${i} \
    --align center \
    --m Final_motifs.meme \
    --p ACRs/ACRs_${i}.fa \
    --n ACRs/ACRs_bkg.fa

done

tomtom \
  -oc tomtom \
  -internal \
  FinalMotifs.meme \
  PlantTFDB_Ath_TF_binding_motifs.meme


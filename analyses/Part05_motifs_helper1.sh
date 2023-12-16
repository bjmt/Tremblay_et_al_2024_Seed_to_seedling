# Benjamin J.M. Tremblay, June 2023
#
# Motif analyses

set -ex

for i in C1 C2 C3 C4 C5 C6 ; do

  streme \
    --oc Promoters_streme/${i} \
    --maxw 12 \
    --seed 1234 \
    --align right \
    --n Promoters/Promoters_bkg.fa \
    --p Promoters/Promoters_${i}.fa

done

for i in A1 A2 A3 A4 A5 ; do

  streme \
    --oc ACRs_streme/${i} \
    --maxw 12 \
    --seed 1234 \
    --align center \
    --n ACRs/ACRs_bkg.fa \
    --p ACRs/ACRs_${i}.fa

done


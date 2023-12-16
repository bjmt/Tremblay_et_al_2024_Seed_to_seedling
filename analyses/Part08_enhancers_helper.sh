# Benjamin J.M. Tremblay, June 2023
#
# Calculate enhancer activity

set -ex

GROUPS=("DS" "S24" "S72" "L6" "L26" "L57")

mkdir -p Enhancer_activities

for i in ${GROUPS[@]} ; do

  bigWigAverageOverBed \
    csrnaseq_bw/${i}.merged.plus.bw \
    Final_enhancers_500bp.bed \
    Enhancer_activities/${i}_plus.txt

  bigWigAverageOverBed \
    csrnaseq_bw/${i}.merged.minus.bw \
    Final_enhancers_500bp.bed \
    Enhancer_activities/${i}_minus.txt

done


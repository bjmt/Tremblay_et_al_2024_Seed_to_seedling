# Benjamin J.M. Tremblay, June 2023
#
# Calculate enhancer activity

set -ex

GROUP_NAMES=("DS" "S24" "S72" "L6" "L26" "L57")

mkdir -p Enhancer_activities

for i in ${GROUP_NAMES[@]} ; do

  bigWigAverageOverBed \
    ${i}.merged.csRNA.plus.bw \
    Final_enhancers_500bp.bed \
    Enhancer_activities/${i}_plus.txt

  bigWigAverageOverBed \
    ${i}.merged.csRNA.minus.bw \
    Final_enhancers_500bp.bed \
    Enhancer_activities/${i}_minus.txt

done


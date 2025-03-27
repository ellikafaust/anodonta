#!/bin/bash
#SBATCH --job-name=ld_decay
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1:00:00
#SBATCH --output=logs/%x.o%A
#SBATCH --error=logs/%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch


###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$

out="ld"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink

BED_FILE="anatina"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.25 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --make-bed --out ${out}/${BED_FILE}_maf_thin
awk '{print $7, $5-$2}' ${out}/${BED_FILE}_maf_thin.ld > ${out}/${BED_FILE}_maf_thin_ld_with_distance.txt 

BED_FILE="anatina_N_mac5"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.25 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --make-bed --out ${out}/${BED_FILE}_maf_thin
awk '{print $7, $5-$2}' ${out}/${BED_FILE}_maf_thin.ld > ${out}/${BED_FILE}_maf_thin_ld_with_distance.txt 

BED_FILE="anatina_TI_mac2"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.25 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --make-bed --out ${out}/${BED_FILE}_maf_thin
awk '{print $7, $5-$2}' ${out}/${BED_FILE}_maf_thin.ld > ${out}/${BED_FILE}_maf_thin_ld_with_distance.txt 

#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$

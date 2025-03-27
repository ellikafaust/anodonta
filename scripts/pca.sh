#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00
#SBATCH --output=logs/%x.o%A
#SBATCH --error=logs/%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch

###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$


out="pca"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink

BED_FILE="anatina_ld"
#BED_FILE="anatina_N_mac5_ld"

# create pca and files that allow proper percentage variance calculations
plink --bfile $BED_FILE --allow-extra-chr --pca var-wts --make-rel --out ${out}/$BED_FILE
#plink --bfile $BED_FILE --allow-extra-chr --pca --out $BED_FILE

#Calculate sum of variance from relation covariance matrix:
sum1=$(awk '{sum+=$NR;}END{print sum}' ${out}/${BED_FILE}.rel)

#Calculate percentage variance explained (pve) and write to file
while read line; do
    echo "$(echo "scale=4; $line / $sum1 * 100" | bc)"
done < ${out}/${BED_FILE}.eigenval > ${out}/${BED_FILE}.pve

# calculate missingness to be plotted on PCs
plink --bfile $BED_FILE --allow-extra-chr --missing --out ${out}/$BED_FILE 


#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$



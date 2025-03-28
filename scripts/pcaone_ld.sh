#!/bin/bash
#SBATCH --job-name=ld
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G
#SBATCH --time=24:00:00
#SBATCH --output=%x.o%A
#SBATCH --error=%x.e%A

###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$


source $GDCstack
module load plink


# split vcf file in two for memory usage
BED_FILE="cygnea_exulcerata"


$HOME/software/PCAone/PCAone -b $BED_FILE -k 3 -n 1 --ld-r2 0.1 --ld-bp 100000 --out $BED_FILE

awk '{print $2 }' ${BED_FILE}.ld.prune.out > snp2rm

# Finally filter out SNPs in LD
#  make a bed file with ld filter ( here I also do pca )
plink --bfile $BED_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --exclude snp2rm --make-bed --pca --out ${BED_FILE}_ld

# make vcf file from filtered bed file
plink --bfile --out ${BED_FILE}_ld --allow-extra-chr --recode vcf-iid --out ${BED_FILE}_ld




#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$



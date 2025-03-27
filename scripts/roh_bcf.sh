#!/bin/bash
#SBATCH --job-name=roh_bcf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00
#SBATCH --output=%x.o%A_%a
#SBATCH --error=%x.e%A_%a
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch

###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$

source /cluster/project/gdc/shared/stack/GDCstack.sh

module load bcftools


for VCF_FILE in all_species anatina_N anatina_TI cygnea exulcerata; do
	echo "working with" $VCF_FILE
    VCF_FILE_ALL="all_species_top20.recode.vcf.gz"
    #vcftools --gzvcf $VCF_FILE_ALL --keep ${VCF_FILE}_samples --recode --out ${VCF_FILE}_top20
    #bgzip ${VCF_FILE}_top20.recode.vcf
    #tabix -p vcf ${VCF_FILE}_top20.recode.vcf.gz
    #bcftools +fill-tags ${VCF_FILE}_top20.recode.vcf.gz -Oz -o ${VCF_FILE}_top20_info.vcf.gz -- -t all
    bcftools roh --GTs-only 30 ${VCF_FILE}_top20_info.vcf.gz -O r -o ${VCF_FILE}_top20_roh.txt 
done



#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$



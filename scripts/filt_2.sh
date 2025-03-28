#!/bin/bash
#SBATCH --job-name=filt_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00
#SBATCH --output=logs/%x.o%A
#SBATCH --error=logs/%x.e%A


##################################################################################################################################################################################################
echo "starting ${SLURM_JOB_ID} at $(date)"
##################################################################################################################################################################################################

VCF_FILE="ac_ae_raw_miss5_Q30_DP2_bi_imiss5.vcf.gz"
echo "Processing file: $VCF_FILE"

source $GDCstack
module load vcftools/0.1.16-tc6l6nq
module load bcftools


# filter site based on missingess and mean site depth between 0.05 and 0.95 quantiles
vcftools --gzvcf ${VCF_FILE} --max-missing 0.8 --min-meanDP 3 --max-meanDP 8.86667 --mac 2 --recode-INFO-all --recode --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss8_meanDP_mac2.vcf.gz} 

# filter for stranbias
bcftools filter -e 'INFO/PV4[0]< 0.01' ${VCF_FILE/.vcf.gz/_miss8_meanDP_mac2.vcf.gz} -Oz -o ${VCF_FILE/.vcf.gz/_miss8_meanDP_mac2_SB.vcf.gz}

# Wait for all background processes to complete
wait

echo "All vcftools analyses are done!"

##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################
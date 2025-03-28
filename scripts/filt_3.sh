#!/bin/bash
#SBATCH --job-name=filt_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x.o%A
#SBATCH --error=logs/%x.e%A



##################################################################################################################################################################################################
echo "starting ${SLURM_JOB_ID} at $(date)"
##################################################################################################################################################################################################

VCF_FILE="ac_ae_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB.vcf.gz"

echo "Processing file: $VCF_FILE"


source $GDCstack
module load vcftools/0.1.16-tc6l6nq
module load bcftools



# filter for missingness
vcftools --gzvcf ${VCF_FILE} --max-missing 0.85 --recode-INFO-all --recode --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss85.vcf.gz} 

VCF_FILE=${VCF_FILE/.vcf.gz/_miss85.vcf.gz} 

echo "Finished filtering and calculating individual stats $(date)"
vcftools --gzvcf ${VCF_FILE} --missing-indv --out ${VCF_FILE/.vcf.gz/} 

tail +2 ${VCF_FILE/.vcf.gz/}.imiss | awk '$5>0.3' | cut -f 1 > samples_imiss3

vcftools --gzvcf ${VCF_FILE} --mac 5 --remove samples_imiss3 --recode --recode-INFO-all --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_imiss3_mac5.vcf.gz} 


echo 'Number of samples to be removed'
wc -l samples_imiss3


# Wait for all background processes to complete
wait

echo "All vcftools analyses are done!"

##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################
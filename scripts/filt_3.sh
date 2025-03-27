#!/bin/bash
#SBATCH --job-name=filt_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x.o%A
#SBATCH --error=logs/%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch


##################################################################################################################################################################################################
echo "starting ${SLURM_JOB_ID} at $(date)"
##################################################################################################################################################################################################

# file need to be provided when submitting the job
# eg: sbatch my_script.sh my_input.vcf.gz

# Check if a file is provided
if [ -z "$1" ]; then
    echo "Usage: sbatch my_script.sh <input_file>"
    exit 1
fi

INPUT_FILE=$1

echo "Processing file: $INPUT_FILE"


source $GDCstack
module load vcftools/0.1.16-tc6l6nq
module load bcftools

VCF_FILE=$INPUT_FILE

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
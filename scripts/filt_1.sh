#!/bin/bash
#SBATCH --job-name=filt_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4G
#SBATCH --time=40:00:00
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

# filter for site missingness, mapping quality, minimum depth, biallelic sites and remove samples with missingness above cutoff
vcftools --gzvcf $VCF_FILE --max-missing 0.5 --minQ 30 --minDP 2 --max-alleles 2 --min-alleles 2 --remove samples_imiss5 --recode --recode-INFO-all --stdout | bgzip -c > ${VCF_FILE/.vcf.gz/_miss5_Q30_DP2_bi_imiss5.vcf.gz}

# stats
VCF_FILE2=${VCF_FILE/.vcf.gz/_miss5_Q30_DP2_bi_imiss5.vcf.gz}

# mean coverage across the sites.
vcftools --gzvcf ${VCF_FILE2} --site-mean-depth --out ${VCF_FILE2/.vcf.gz/} &

# missingnes per site
vcftools --gzvcf ${VCF_FILE2} --missing-site --out ${VCF_FILE2/.vcf.gz/} &

# Wait for all background processes to complete
wait

echo "All vcftools analyses are done!"

##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################
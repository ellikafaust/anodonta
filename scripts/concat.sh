#!/bin/bash
#SBATCH --job-name=concat
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

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load bcftools
module load vcftools/0.1.16-tc6l6nq

VCF_FILE="ac_ae_raw.vcf.gz"

# make list of bcfs to combine
ls ./SNPs/*.bcf > bcf_list

# Concatenate bcf files
bcftools concat --file-list bcf_list -Oz -o $VCF_FILE

echo "Finished combining. Start calculating stats $(date)"

# missingness per individual
vcftools --gzvcf ${VCF_FILE} --missing-indv --out ${VCF_FILE/.vcf.gz/} &

# mean coverage across the sites.
vcftools --gzvcf $VCF_FILE --site-mean-depth --out ${VCF_FILE/.vcf.gz/} &

# individual depth
vcftools --gzvcf ${VCF_FILE} --depth --out ${VCF_FILE/.vcf.gz/} &

# Wait for all background processes to complete
wait

echo "All vcftools analyses are done!"

##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################
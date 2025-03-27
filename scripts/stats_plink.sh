#!/bin/bash
#SBATCH --job-name=stats
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
#SBATCH --output=logs/stats.o%A
#SBATCH --error=logs/stats.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch


###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$

module load stack
module load plink

BED_FILE="all_species"


plink --bfile $BED_FILE --allow-extra-chr --het --out $BED_FILE &
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_samples --recode vcf-iid  --out anatina

echo "Finished filtering and calculating individual stats $(date)"
vcftools --gzvcf $VCF_FILE --depth --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --missing-indv --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --het --out ${VCF_FILE/.vcf.gz/}
#vcftools --gzvcf $VCF_FILE --relatedness --out ${VCF_FILE/.vcf.gz/}
#vcftools --gzvcf $VCF_FILE --relatedness2 --out ${VCF_FILE/.vcf.gz/}

echo "Calculating allele stats $(date)"
vcftools --gzvcf $VCF_FILE --freq2 --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --counts2 --out ${VCF_FILE/.vcf.gz/}

echo "Calculating site stats $(date)"
#vcftools --gzvcf $VCF_FILE --site-depth --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --site-mean-depth --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --missing-site --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --singletons --out ${VCF_FILE/.vcf.gz/}
#vcftools --gzvcf $VCF_FILE --hardy --out ${VCF_FILE/.vcf.gz/}
##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################


#!/bin/bash
#SBATCH --job-name=stats
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/stats.o%A
#SBATCH --error=logs/stats.e%A
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
OUT_PREFIX="${VCF_FILE/.vcf.gz/}"
OUT_DIR="stats_${OUT_PREFIX}"

# if there is no out directory make one
if [ ! -e ${OUT_DIR} ]  ; then mkdir ${OUT_DIR} ; fi

echo "Finished filtering and calculating individual stats $(date)"
vcftools --gzvcf $VCF_FILE --depth --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --missing-indv --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --het --out $OUT_DIR/$OUT_PREFIX &
#vcftools --gzvcf $VCF_FILE --relatedness --out $OUT_DIR/$OUT_PREFIX &
#vcftools --gzvcf $VCF_FILE --relatedness2 --out $OUT_DIR/$OUT_PREFIX &

echo "Calculating allele stats $(date)"
vcftools --gzvcf $VCF_FILE --counts2 --out $OUT_DIR/$OUT_PREFIX &

echo "Calculating site stats $(date)"
#vcftools --gzvcf $VCF_FILE --site-depth --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --site-mean-depth --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --missing-site --out $OUT_DIR/$OUT_PREFIX &
#vcftools --gzvcf $VCF_FILE --singletons --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --hardy --out $OUT_DIR/$OUT_PREFIX &

wait
##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################


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


VCF_FILE=$INPUT_FILE
BED_FILE="cygnea_exulcerata"
OUT_DIR="beds"

# if there is no out directory make one
if [ ! -e ${OUT_DIR} ]  ; then mkdir ${OUT_DIR} ; fi

source $GDCstack
module load vcftools/0.1.16-tc6l6nq
module load bcftools
module load plink

#  make a bed file
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out ${OUT_DIR}/${BED_FILE}

# make a vcf file which has less info field and thus is smaller and easier to read into memory
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --recode vcf-iid --out ${OUT_DIR}/$BED_FILE

###### CYGNEA #########

awk -F',' '$1 == "Anodonta cygnea" {print $NF "\t" $NF}' metadata_all.csv > cygnea_samples
plink --bfile $BED_FILE --allow-extra-chr --keep cygnea_samples --make-bed --pca --out ${OUT_DIR}/cygnea

grep -vE "AC83|AC90" cygnea_samples > cygnea_samples_exc_hy
plink --bfile $BED_FILE --allow-extra-chr --keep cygnea_samples_exc_hy --mac 5 --make-bed --pca --recode vcf-iid --out ${OUT_DIR}/cygnea_mac5 


###### EXULCERATA #########

awk -F',' '$1 == "Anodonta exulcerata" {print $NF "\t" $NF}' metadata_all.csv  > exulcerata_samples
plink --bfile $BED_FILE --allow-extra-chr --keep exulcerata_samples --make-bed --pca --out ${OUT_DIR}/exulcerata 
plink --bfile $BED_FILE --allow-extra-chr --keep exulcerata_samples --mac 2 --make-bed --pca --recode vcf-iid  --out ${OUT_DIR}/exulcerata_mac2 


wait
##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################

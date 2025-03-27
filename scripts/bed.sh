#!/bin/bash
#SBATCH --job-name=bed
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x.o%A
#SBATCH --error=logs/%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch


###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$

module load stack
module load plink
module load vcftools/0.1.16-kpfcvxs
module load bcftools

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac5_imiss30.vcf.gz"
BED_FILE="all_species"

#  make a bed file
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out ${BED_FILE}

# make a vcf file which has less info field and thus is smaller and easier to read into memory
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --recode vcf-iid --out $BED_FILE

# to create a tab delimited file for plink of individuals to keep
awk -F',' '$1 == "Anodonta anatina" {print $NF "\t" $NF}' groups_aa.csv  > anatina_samples
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_samples --make-bed --pca --out anatina
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_samples --recode vcf-iid  --out anatina


###### CYGNEA #########

awk -F',' '$1 == "Anodonta cygnea" {print $NF "\t" $NF}' groups_aa.csv  > cygnea_samples
plink --bfile $BED_FILE --allow-extra-chr --keep cygnea_samples --make-bed --pca --out cygnea
plink --bfile $BED_FILE --allow-extra-chr --keep cygnea_samples --recode vcf-iid  --out cygnea

###### EXULCERATA #########

awk -F',' '$1 == "Anodonta exulcerata" {print $NF "\t" $NF}' groups_aa.csv  > exulcerata_samples
plink --bfile $BED_FILE --allow-extra-chr --keep exulcerata_samples --make-bed --pca --out exulcerata
plink --bfile $BED_FILE --allow-extra-chr --keep exulcerata_samples --recode vcf-iid  --out exulcerata

###### EXULCERATA & CYGNEA #########
plink --bfile $BED_FILE --allow-extra-chr --remove anatina_samples --make-bed --pca --out exulcerata_cygnea
plink --bfile $BED_FILE --allow-extra-chr --remove anatina_samples --recode vcf-iid  --out exulcerata_cygnea


#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$

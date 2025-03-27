#!/bin/bash
#SBATCH --job-name=hwe
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=%x.o%A
#SBATCH --error=%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch


###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load bcftools
module load  vcftools/0.1.16-tc6l6nq

# Individual heterozygosity to calculate FIS
VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac5_imiss30.vcf.gz"

awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 {print $NF "\t" $NF "\t" $2 "_" $5}' groups_aa.csv > pheno

# get list of individuals which passed filtering
bcftools query -l $VCF_FILE > individuals
# Extract unique groups (populations) from the third column of the pheno file
grep -Ff individuals pheno | cut -f 3 | sort | uniq > pops



while read -r group; do
    # Create the 'keep' file for the group
    grep -Ff individuals pheno | awk -v grp="$group" '$3 == grp {print $1}' > keep$group
    
    # Check the number of lines (individuals) in the 'keep' file
    num_individuals=$(wc -l < keep$group)

    # Only run vcftools if there are more than 5 individuals
    if [ "$num_individuals" -gt 9 ]; then
        echo "Running vcftools for group $group with $num_individuals individuals"
        vcftools --gzvcf $VCF_FILE --keep keep$group --hardy --out $group
    else
        echo "Skipping group $group, not enough individuals ($num_individuals)"
    fi
done < pops


#awk -F',' '$1 == "Anodonta cygnea" {print $NF}' groups_aa.csv | grep -vE "AC83|AC90" -- | grep -Ff individuals -- > cygnea_samples_exc_hy
#awk -F',' '$1 == "Anodonta anatina" && $6 != "Ticino" {print $NF}' groups_aa.csv | grep -Ff individuals > anatina_N_samples
#vcftools --gzvcf $VCF_FILE --keep anatina_N_samples --hardy --out anatina_N &
#vcftools --gzvcf $VCF_FILE --keep cygnea_samples_exc_hy --hardy --out cygnea &

#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$

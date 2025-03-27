#!/bin/bash
#SBATCH --job-name=info
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
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


VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac5_imiss30.vcf.gz"


awk -F',' '$1 == "Anodonta cygnea {print $NF}' groups_aa.csv  > cygnea_samples
#awk -F',' '$1 == "Anodonta cygnea" && $5 == "UFS" {print $NF}' groups_aa.csv  > cygnea_UFS

bcftools view -S cygnea_samples --force-samples -o ${VCF_FILE/.vcf.gz/_cygnea.vcf} -Oz ${VCF_FILE}

# tabix -p vcf ${VCF_FILE/.vcf.gz/_cygnea.vcf.gz}

bcftools +fill-tags  -o ${VCF_FILE/.vcf.gz/_ac.vcf.gz} -Oz ${VCF_FILE/.vcf.gz/_cygnea.vcf.gz} -- -t all

#bcftools index ${VCF_FILE/.vcf.gz/_ac.vcf.gz}


bcftools query -f '%CHROM\t%POS\t%INFO/MQ\t%INFO/DP\t%INFO/FS\n' AC_top20_phased.bcf | head

bcftools query -f '%CHROM\t%POS\t%INFO/MQ\t%INFO/DP\t%INFO/FS\n' ${VCF_FILE}_cygnea_updated | head


awk 'BEGIN {OFS="\t"} {chrom=$1; mq=$2; fs=$3; dp=$4; 
sum_mq[chrom]+=mq; sum_fs[chrom]+=fs; sum_dp[chrom]+=dp; count[chrom]++} 
END {for (chrom in sum_mq) print chrom, sum_mq[chrom]/count[chrom], sum_fs[chrom]/count[chrom], sum_dp[chrom]/count[chrom]}' variant_info.txt > mean_per_chromosome.txt


awk -F',' '$1 == "Anodonta cygnea" {print $NF}' groups_aa.csv  > cygnea_samples

bcftools view -S anatina_samples -o ${VCF_FILE}_anatina $VCF_FILE

bcftools view -S cygnea_samples $VCF_FILE | bcftools +fill-tags  -o ${VCF_FILE}_cygnea -O v -- -t MQ,FS,DP


bcftools view -h ${VCF_FILE}_cygnea | grep "##INFO"


bcftools query -f '%CHROM\t%POS\t%INFO/MQ\t%INFO/DP\t%INFO/FS\n' updated.vcf | head


###############################################################################################################################################$
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$


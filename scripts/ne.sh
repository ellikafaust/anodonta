#!/bin/bash
#SBATCH --job-name=ne
#SBATCH --array=1%1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=500
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x.o%A_%a
#SBATCH --error=logs/%x.e%A_%a


############################################################################################################################$
echo "starting ${SLURM_JOB_ID} at $(date)"
############################################################################################################################$


VCF_FILE="exulcerata_mac2"
pops=pops_10_AE
out=ne


# how to reference each sample by array number
IDX=$SLURM_ARRAY_TASK_ID
group=`sed -n ${IDX}p < $pops`

echo "Processing file: $pops"

source $GDCstack
module load bcftools

# Create the 'keep' file for the group
#awk -v grp="$group" '$3 == grp {print $1}' pheno > ${out}/keep$group
    
#bcftools view -T ${out}/${VCF_FILE}_2Msnsp.txt -S ${out}/keep${group} -o ${out}/${group}.vcf -O v ${VCF_FILE}.vcf.gz

echo "Staring Ne calculations of $pops"

$HOME/software/currentNe/currentNe -t 8 ${out}/${group}.vcf 19



##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################
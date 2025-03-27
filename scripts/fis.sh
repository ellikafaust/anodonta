#!/bin/bash
#SBATCH --job-name=fis
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

source $GDCstack
module load vcftools/0.1.16-tc6l6nq

out="het"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

BED_FILE="anatina"

# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 && $9 == "passed" {print $NF "\t" $NF "\t" $2 "_" $5}' metadata_all.csv > pheno
# Extract unique groups (populations) from the third column of the pheno file
tail -n +2 pheno | cut -f 3 | sort | uniq  > pops

grep "AA_" pops > pops_AA
# Individual heterozygosity to calculate FIS


while read -r group; do
    echo "working with group $group"
    awk -v grp="$group" '$3 == grp {print $1 "\t" $2 }' pheno > ${out}/keep$group
    #plink --bfile $BED_FILE --allow-extra-chr --keep keep$group --het --out ${BED_FILE}${group}
    vcftools --vcf ${BED_FILE}.vcf --keep ${out}/keep$group --het --out ${out}/${BED_FILE}${group}
done < pops_AA

# concatenate the files


echo -e "pop\tsample\tO(HOM)\tE(HOM)\tN_SITES\tF" > ${out}/${BED_FILE}_pop_het.txt

while read -r group; do
tail -n +2 ${out}/*${group}.het | awk -v grp="$group" '{print grp "\t" $0}' >> ${out}/${BED_FILE}_pop_het.txt
done < pops_AA



#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$

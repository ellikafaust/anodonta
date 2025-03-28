#!/bin/bash
#SBATCH --job-name=fis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x.o%A
#SBATCH --error=logs/%x.e%A


###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$

source $GDCstack
module load vcftools/0.1.16-tc6l6nq

out="het"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

BED_FILE="cygnea_exulcerata"



# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 && $9 == "passed" {print $NF "\t" $NF "\t" $2 "_" $5}' metadata_all.csv > pheno
# Extract unique groups (populations) from the third column of the pheno file
tail -n +2 pheno | cut -f 3 | sort | uniq | grep -v "AA_" > pops

grep -v "AA_" pops > pops_AC_AE
# Individual heterozygosity to calculate FIS


while read -r group; do
    echo "working with group $group"
    awk -v grp="$group" '$3 == grp {print $1 "\t" $2 }' pheno > ${out}/keep$group
    #plink --bfile $BED_FILE --allow-extra-chr --keep keep$group --het --out ${BED_FILE}${group}
    vcftools --vcf ${BED_FILE}.vcf --keep ${out}/keep$group --het --out ${out}/${BED_FILE}${group}
done < pops_AC_AE

# concatenate the files


echo -e "pop\tsample\tO(HOM)\tE(HOM)\tN_SITES\tF" > ${BED_FILE}_pop_het.txt

while read -r group; do
tail -n +2 *${group}.het | awk -v grp="$group" '{print grp "\t" $0}' >> ${BED_FILE}_pop_het.txt
done < pops_AC_AE




#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$

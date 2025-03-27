#!/bin/bash
#SBATCH --job-name=pi
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

out="pi"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load vcftools/0.1.16-tc6l6nq

BED_FILE="anatina"

# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 && $9 == "passed" {print $NF "\t" $NF "\t" $2 "_" $5}' metadata_all.csv > pheno
# Extract unique groups (populations) from the third column of the pheno file
tail -n +2 pheno | cut -f 3 | sort | uniq > pops

grep "AA_" pops > pops_AA


while read -r group; do
    # Create the 'keep' file for the group
    awk -v grp="$group" '$3 == grp {print $1}' pheno > ${out}/keep$group
    
    # Check the number of lines (individuals) in the 'keep' file
    num_individuals=$(wc -l < ${out}/keep$group)

    # Only run vcftools if there are more than 5 individuals
    if [ "$num_individuals" -gt 4 ]; then
        echo "Running vcftools for group $group with $num_individuals individuals"
        vcftools --vcf ${BED_FILE}.vcf --window-pi 10000 --keep ${out}/keep$group --out ${out}/pi_10kb_$group
    else
        echo "Skipping group $group, not enough individuals ($num_individuals)"
        rm ${out}/keep$group
    fi
done < pops_AA


# concatenate the files
echo -e "pop\tCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI" > ${out}/${BED_FILE}_pi_10kb.txt
while read -r group; do
tail -n +2 ${out}/pi_10kb_${group}.windowed.pi | awk -v grp="$group" '{print grp "\t" $0}' >> ${out}/${BED_FILE}_pi_10kb.txt
done < pops_AA


#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$

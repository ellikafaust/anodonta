#!/bin/bash
#SBATCH --job-name=dxy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=4:00:00
#SBATCH --output=%x.o%A
#SBATCH --error=%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch


###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load python

echo -e "Pop1\tPop2\tAverage_dxy" > all_dxy_results.txt

# Pairwise comparison and appending results
for file1 in *.afreq; do
    for file2 in *.afreq; do
        if [[ "$file1" < "$file2" ]]; then
            echo "Calculating dxy for $file1 and $file2"
            # Capture the output of the Python script
            dxy=$(python calculate_average_dxy.py "$file1" "$file2" | grep "Average dxy" | awk '{print $3}')
            # Append to the results file
            echo -e "${file1%.afreq}\t${file2%.afreq}\t$dxy" >> all_dxy_results.txt
        fi
    done
done

SNP=26704705
REF=2896208859
FRACTION=$(echo "scale=10; $SNP / $REF" | bc)
# .0092205729

awk -v frac="$FRACTION" '{print $0, $NF * frac}' all_dxy_results.txt > updated_all_dxy_results.txt


#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$

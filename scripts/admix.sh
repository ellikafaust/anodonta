#!/bin/bash
#SBATCH --job-name=admix
#SBATCH --array=1-12%3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --time=10:00:00
#SBATCH --output=admix.o%A_%a
#SBATCH --error=admix.e%A_%a
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch

echo "starting ${SLURM_JOB_ID} at $(date)"
##################################################################################################################################################################################################

# file need to be provided when submitting the job
# eg: sbatch my_script.sh my_input

# Check if a file is provided
if [ -z "$1" ]; then
    echo "Usage: sbatch my_script.sh <input_file>"
    exit 1
fi

INPUT_FILE=$1

echo "Processing file: $INPUT_FILE"

BED_FILE=$INPUT_FILE

source $GDCstack
module load admixture

# how to reference each sample by array number
K=$SLURM_ARRAY_TASK_ID

admixture --cv=5 -j8 ${BED_FILE}.bed $K | tee log_${BED_FILE}_${K}.out

#################################
echo "Job: K=${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$



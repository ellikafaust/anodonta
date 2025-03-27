#!/bin/bash
#SBATCH --job-name=bcall_loop
#SBATCH --array=1-20%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x.o%A_%a
#SBATCH --error=logs/%x.e%A_%a
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch


##################################################################################################################################################################################################
echo "starting ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} at $(date)"
##################################################################################################################################################################################################


source /cluster/project/gdc/shared/stack/GDCstack.sh
module load bcftools

out=SNPs
ref=ref/v3.asm.bp.p_ctg.fa
bams=bam_list
chunks=scaffold_chunks/scaffold_chunk_

# if there is no out directory make one
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

#down-sample region with > max. cov. per individual 
maxcov=50

IDX=$SLURM_ARRAY_TASK_ID

# Set the maximum number of concurrent jobs (adjust based on your CPU cores)
max_jobs=4  # Number of concurrent jobs (adjust as needed)
current_jobs=0  # Counter for active jobs

# Loop through each region/sample in the input file
while IFS= read -r name; do
    echo "Starting calling $name at $(date)"

    # Run the command in the background to process multiple regions concurrently
    bcftools mpileup -f ${ref} --skip-indels -r ${name} -d ${maxcov} -a 'FORMAT/DP,FORMAT/AD,INFO/SCR' -b ${bams}| \
    bcftools call -mv -a GQ,PV4 -Ob -o ${out}/raw.${name}.bcf &

    ((current_jobs++))

    # If we reach the maximum number of concurrent jobs, wait for one to finish
    if ((current_jobs >= max_jobs)); then
        wait -n  # Wait for any job to finish
        ((current_jobs--))
    fi

done < ${chunks}${IDX}

# Wait for all remaining jobs to finish before exiting the script
wait

##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################
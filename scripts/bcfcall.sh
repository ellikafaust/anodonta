#!/bin/bash
#SBATCH --job-name=bcall
#SBATCH --array=1-100%20
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=1G
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
chunks=scaffold_chunks/largest_scaffolds_100

# if there is no out directory make one
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

IDX=$SLURM_ARRAY_TASK_ID
name=`sed -n ${IDX}p < $chunks`


#down-sample region with > max. cov. per individual 
maxcov=50

bcftools mpileup -f ${ref} --skip-indels -r ${name} -d ${maxcov} -a 'FORMAT/DP,FORMAT/AD,INFO/SCR' -b ${bams} | bcftools call -mv -a GQ,PV4 -Ob  -o ${out}/raw.${name}.bcf
# bcftools mpileup -f ${Ref} --skip-indels -r ${name} -d ${maxcov} -a AD,DP,SP,SCR -b $bams -Ou | bcftools call -mv -a GQ,GP,PV4 -Ob  -o ${out}/raw.${name}.bcf

##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################
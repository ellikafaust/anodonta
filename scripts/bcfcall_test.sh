#!/bin/bash
#SBATCH --job-name=bcall
#SBATCH --array=1-3%3
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=1G
#SBATCH --time=24:00:00
#SBATCH --output=bcall.o%A_%a
#SBATCH --error=bcall.e%A_%a
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch




##################################################################################################################################################################################################
echo "starting ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} at $(date)"
##################################################################################################################################################################################################

#copy bam files to scratch and generate list
#ls $PANA/bams/*bam > bam_list
#awk -F',' '$1 == "Anodonta anatina" {print $NF "_sort_Q20_nodup.bam"}' groups_aa.csv  | grep -Ff - bam_list > anatina_bams
#awk -F',' '$1 == "Anodonta cygnea" {print $NF "_sort_Q20_nodup.bam"}' groups_aa.csv | grep -Ff - bam_list > cygnea_bams
#awk -F',' '$1 == "Anodonta exulcerata" {print $NF "_sort_Q20_nodup.bam"}' groups_aa.csv | grep -Ff - bam_list > exulcerata_bams
#ls ./*_bams > ind_list

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load bcftools

out=/cluster/scratch/efaust/gt_calling/SNPs
Ref=/cluster/scratch/efaust/gt_calling/ref/Aanat_demux_hifi_l3.bp.p_ctg.fa
pos="common_bad_pos.txt"

# if there is no out directory make one
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

IDX=$SLURM_ARRAY_TASK_ID
name=`sed -n ${IDX}p < ind_list`


#down-sample region with > max. cov. per individual 
maxcov=30

bcftools mpileup -f ${Ref} --skip-indels -R ${pos} -d ${maxcov} -a 'FORMAT/DP,FORMAT/AD' -b ${name} | bcftools call -mv -Ob -o ${name}_raw.bcf
# bcftools mpileup -f ${Ref} --skip-indels -r ${name} -d ${maxcov} -a AD,DP,SP,SCR -b $bams -Ou | bcftools call -mv -a GQ,GP,PV4 -Ob  -o ${out}/raw.${name}.bcf

##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################

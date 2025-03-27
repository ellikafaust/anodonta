#!/bin/bash
#SBATCH --job-name=rm_dup  
#SBATCH --array=1-168%11
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G  
#SBATCH --time=12:00:00     
#SBATCH --tmp=40G
#SBATCH --output=logs/%x.o%A_%a
#SBATCH --error=logs/%x.e%A_%a
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch

##################################################################################################################################################################################################
echo "Starting ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} at $(date)"
##################################################################################################################################################################################################

# if you didn't remove pcr duplicates as part of the mapping job.

#####Definition of the paths:
data=${PANA}/fastqs # where the fastq files are kept
ref=${SCRATCH}/mapping/ref/v3.asm.bp.p_ctg.fa.gz # the reference genome
out=${SCRATCH}/mapping/mapping_stats # where to put the mapping statistics
proj=${PANA}/bams_cygnea # final destination of the bams
bams=${SCRATCH}/bams_cygnea # temporary location for bams in the scratch for later use for snp-calling
sample_list=${SCRATCH}/mapping/ac_ae_sample_list

##Number_processors per Job, hyper-threading activated request 2 but provide BWA 4 CPUs
cpu=4
mem=5GB

##Mapping Quality
Qual=20

# use old softwares tack
source $GDCstack
module load gcc 
module load bwa-mem2
module load sambamba
module load picard

# how to reference each sample by array number
IDX=$SLURM_ARRAY_TASK_ID
name=`sed -n ${IDX}p < ${sample_list}`

# create directories
if [ ! -e ${out}/statsQ${Qual}_dup ]  ; then mkdir ${out}/statsQ${Qual}_dup ; fi
if [ ! -e ${out}/stats_dup ]  ; then mkdir ${out}/stats_dup ; fi

# 
echo "Start processing ${name}"

# mark and remove pcr duplicates
echo "Remove PCR duplicates"

# Some cluster environments provide a picard wrapper script instead of using Java directly. 
# Unfortunately, this wrapper does not accept -Xmx options directly. 
#If your system uses this, you may need to set the memory limit with an environment variable before running Picard:
export _JAVA_OPTIONS="-Xmx5G"

picard MarkDuplicates \
	-TMP_DIR ${TMPDIR} \
	-INPUT ${bams}/${name}_sort_Q${Qual}.bam \
	-OUTPUT ${bams}/${name}_sort_Q${Qual}_nodup.bam \
	-METRICS_FILE ${out}/stats_dup/${name}_nodup \
	-VALIDATION_STRINGENCY LENIENT \
	-MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 512 \
	-REMOVE_DUPLICATES true

# remove bam pre-deduplication
rm ${bams}/${name}_sort_Q${Qual}.bam*

# index bam
sambamba index -t ${cpu} ${bams}/${name}_sort_Q${Qual}_nodup.bam

# mapping statistics after removing pcr duplicates
sambamba flagstat -t ${cpu} ${bams}/${name}_sort_Q${Qual}_nodup.bam > ${out}/statsQ${Qual}_dup/${name}nodup_Q20_flagstat

# copy bam and bai to project folder and remove from SCRATCH
cp ${bams}/${name}_sort_Q${Qual}_nodup.bam* ${proj}/

# remove from scratch (if you are not doing snp-calling straight away)
#rm ${bams}/${name}_sort_Q${Qual}.bam*

##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################


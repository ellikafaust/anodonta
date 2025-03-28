#!/bin/bash
#SBATCH --job-name=map 
#SBATCH --array=1-168%11
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G  
#SBATCH --time=24:00:00     
#SBATCH --tmp=40G
#SBATCH --output=%x.o%A_%a
#SBATCH --error=%x.e%A_%a

##################################################################################################################################################################################################
echo "Starting ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} at $(date)"
##################################################################################################################################################################################################

# make a list concatenated samples
#ls $PANA/fastqs/*_R1_001.fastq.gz | sed 's/^[^-]*-\([^-]*\)_S.*$/\1/' > sample_list
#389 samples
# only keeping samples which are not Anodonta anatina (ie cygnea and exulcerata)
#awk -F',' '$1 != "Anodonta anatina" {print $NF}' groups_aa.csv | sed 's/^[^-]*-\([^-]*\)_S.*$/\1/' | grep -Ff - sample_list > ac_ae_sample_list


#####Definition of the paths:
data=${PANA}/fastqs # where the fastq files are kept
ref=${SCRATCH}/mapping/ref/Acygnea.v3.asm.bp.p_ctg.fa.gz # the reference genome
out=${SCRATCH}/mapping/mapping_stats # where to put the mapping statistics
proj=${PANA}/bams_cygnea # final destination of the bams
bams=${SCRATCH}/bams_cygnea # temporary location for bams in the scratch for later use for snp-calling
sample_list=${SCRATCH}/mapping/ac_ae_sample_list # sample list

##Number_processors per Job, hyper-threading activated request 2 but provide BWA 4 CPUs
cpu=4
mem=5GB

##Mapping Quality
Qual=20

# load softwares
source $GDCstack
module load gcc 
module load bwa-mem2
module load sambamba
module load picard

# how to reference each sample by array number
IDX=$SLURM_ARRAY_TASK_ID
name=`sed -n ${IDX}p < ${sample_list}`

echo "Start processing ${name}"

###copy data to node tmpdir
cp ${data}/*-${name}_*.fastq.gz ${TMPDIR}/

#re-name files
mv ${TMPDIR}/*-${name}_*R1*.fastq.gz ${TMPDIR}/${name}_R1.fq.gz
mv ${TMPDIR}/*-${name}_*R2*.fastq.gz ${TMPDIR}/${name}_R2.fq.gz

# create folder for mapping statistics
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi
if [ ! -e ${out}/stats ]  ; then mkdir ${out}/stats ; fi
if [ ! -e ${out}/statsQ${Qual} ]  ; then mkdir ${out}/statsQ${Qual} ; fi
if [ ! -e ${out}/stats_dup ]  ; then mkdir ${out}/stats_dup ; fi

# create folders for bams
if [ ! -e ${proj} ]  ; then mkdir ${proj} ; fi
if [ ! -e ${bams} ]  ; then mkdir ${bams} ; fi

# map reads to genome 
echo "Map paired-end read using default parameters"
bwa-mem2 mem ${ref} ${TMPDIR}/${name}_R1.fq.gz ${TMPDIR}/${name}_R2.fq.gz -M -R "@RG\tID:${name}\tSM:${name}\tPL:Illumina" -t ${cpu} > ${TMPDIR}/${name}.sam

#remove fastq file from temp
rm ${TMPDIR}/${name}*.gz

# sort sam and convert to bam
echo "Convert sam to bam and sort it"
sambamba view -t ${cpu} -S ${TMPDIR}/${name}.sam -f bam -o /dev/stdout | sambamba sort /dev/stdin -o ${TMPDIR}/${name}_sort.bam -t ${cpu} -m ${mem} --tmpdir ${TMPDIR}

###delete sam file
rm ${TMPDIR}/${name}.sam

# mapping statistics before any filtering
echo "Mapping statistics"
sambamba flagstat -t ${cpu} ${TMPDIR}/${name}_sort.bam > ${out}/stats/${name}_flagstat

# filtering for mapping quality
echo "Remove low quality reads"
sambamba view  -F "mapping_quality >= "${Qual} ${TMPDIR}/${name}_sort.bam -o ${bams}/${name}_sort_Q${Qual}.bam -t ${cpu} -f bam

# remove unfiltered bams
rm ${TMPDIR}/${name}_sort.bam 

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

# copy bam and bai to project folder
cp ${bams}/${name}_sort_Q${Qual}_nodup.bam* ${proj}/

# remove from scratch (if you are not doing snp-calling straight away)
#rm ${bams}/${name}_sort_Q${Qual}.bam*


##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
echo "Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################


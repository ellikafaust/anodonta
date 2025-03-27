# WGS pipeline for duckmussel Anodonta anatina
# Last updated by Ellika Faust September 2024

# =========================================================================
### 0. Useful commands etc.
# =========================================================================

# user guides
# https://www.gdc-docs.ethz.ch/EulerManual/site/
# https://www.gdc-docs.ethz.ch/GeneticDiversityAnalysis/GDA21/site/snps/

# available submit scripts
/cluster/project/gdc/shared/scripts/submitscripts/


P976=/cluster/work/gdc/shared/p976

/cluster/work/gdc/shared/p976


# submit script 
######### !OBS! ###########
# -n -c and -p is used differently to Rackham!

#SBATCH --ntasks=1               #Requesting 1 node (is always 1)
#SBATCH --cpus-per-task=1        #Requesting 1 CPU
#SBATCH --mem-per-cpu=2G         #Requesting 2 Gb memory per core 
# example
sbatch -n 1 -c 1 --mem-per-cpu 2G -t 1:00:00 -J s2b -o s2b.%J.out -e s2b.e%J.err --mail-type ALL --mail-user julie.conrads@eawag.ch s2b.sh

# interactive
srun --cpus-per-task=1 --time=04:00:00 --pty bash #
srun --cpus-per-task=1 --mem-per-cpu=4G --time=04:00:00 --pty bash #

# Overview of the submitted jobs (pending and running)
squeue

#chain a job
sbatch --dependency=afterok:"46888211_9:46888211_10" < jobB.slurm.sh

# CPU and memory usage of running jobs
myjobs -r
myjobs -j <Job-ID>/<Array-ID>

# Efficiency of finished job(s)
source $GDCstack
module load reportseff/2.7.6
reportseff 18527789
reportseff 18580780
reportseff 18683418


# Connect to a node to check on a job. This is an advanced command. Can be useful to check on the real-time CPU usage or data in the local scratch space.
srun --interactive --jobid <job-ID> --pty bash

# Disk usage of your personal home and your scratch
/cluster/apps/local/lquota

# Available disk space of the GDC share (only works if you are on the works file system)
/cluster/apps/local/lquota /cluster/work/gdc/shared/p976

# Archive tar.gz
tar cvzf <folder>.tar.gz <folder>
# Extract archive
tar xvf <folder>.tar.gz
# list content of archive
tar -ztvf file.tar.gz
# Extract specific file
tar -zxvf file.tar.gz path/to/file

# Archive fa.gz
gzip -d <folder>.fa.gz

# Project directories
/cluster/work/gdc/shared/p979 # Grayling
/cluster/work/gdc/shared/p976 # mussels

# module stacks
lmod2env #new software stacks
env2lmod #old software stack
# source /cluster/apps/local/lmod2env

source $GDCstack
module load vcftools/0.1.16-kpfcvxs



# ==================================================================================================================================================
# Outline
# ==================================================================================================================================================
# 1. Download genome and index
# 2. Download data
# 3. Quality check (FastQC)
# 4. Mapping reads to genome (bwa-mem2)
# 5. Variant calling (bcftools)
# 6. Filtering (vcftools)
# 7. Data analysis



# ==================================================================================================================================================
## 1. Download genome
# ==================================================================================================================================================

# -------------------------------------------------------------------------
# Index genome

# index genome with bwa-mem2
echo '#!/bin/bash
#SBATCH --job-name=indexv4 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=1:00:00
#SBATCH --output=%x.o%A
#SBATCH --error=%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch

source $GDCstack
module load bwa-mem2
bwa-mem2 index ./ref/v4.asm.bp.p_ctg.fa.gz
'> index4.sh

sbatch index4.sh
#Total time taken: 1848.6002

# --------  copy genomes and indexes to project ----------

# ==================================================================================================================================================
# 2. Download data
# ==================================================================================================================================================

# -------------------------------------------------------------------------
# Download data from fgcz gstore with Alex Login into project folder
#Test with single file
wget --user weberalexandraanhthu --ask-password https://fgcz-gstore.uzh.ch/projects/p30854/o32826_NovaSeq_230928_X007/328261_05-AA05_S17_R2_001.fastq.gz

#Download full folder
wget --user weberalexandraanhthu -e robots=off --ask-password -r --no-parent -nH --cut-dirs=2 --reject="index.html*" https://fgcz-gstore.uzh.ch/projects/p30854/o32826_NovaSeq_230928_X007/


# ==================================================================================================================================================
# 3. Quality check
# ==================================================================================================================================================

# checking all the data can take some time so it can be good to copy over a sub-sample to the SCRATCH system and check the quality of them

# -------- copy fastq.gz files to SCRATCH ----------

# make a list concatenated samples
ls $PANA/fastqs/*_R1_001.fastq.gz | sed 's/^[^-]*-\([^-]*\)_S.*$/\1/' > sample_list
#389 samples

# only keeping samples with are not Anodonta anatina (ie cygnea and exulcerata)
awk -F',' '$1 != "Anodonta anatina" {print $NF}' metadata_all.csv | grep -Ff - sample_list > ac_ae_sample_list
# 168 samples (159 cygnea and 17 exulcerata)

# open fastqc script and adjust paths and array length
nano fastqc.sh

sbatch fastqc.sh

# make a summary report of fastqc
module load multiqc
multiqc ./


# --------  tar archive and copy fastqc, multiqc and log files to project ----------

# ==================================================================================================================================================
# 4. Mapping
# ==================================================================================================================================================


# make a sample list and copy
ls *_R1_001.fastq.gz |sed 's/^[^-]*-\([^-]*\)_S.*$/\1/'> sample_list


cp ${data}/*.fastq.gz

# Open script and change paths and names
nano map.sh

# map.sh will do the following for one sample at a time (in arrays)
# 1. concatenate lanes
# 2. map reads with bcftools
# 3. sort and convert sam to bam
# 4. calculate mapping statistics before filtering
# 5. filter for Q20 and remove pcr duplicates
# 6. index and calculate statistics after filtering
# 7. copy filtered bams to the project and delete all intermediate files
sbatch map.sh

# after completed run check all log files for potential error messages and samples which did not complete

# summaries flagstat and pcr de-duplication log files using multiqc
source $GDCstack
module load multiqc

multiqc mapping_stats/stats/ -n multiqc
multiqc mapping_stats/stats_dup/ -n multiqc_dup
multiqc mapping_stats/statsQ20/ -n multiqc_Q20
multiqc mapping_stats/statsQ20_dup/ -n multiqc_Q20_dup


# ==================================================================================================================================================
# 6. Filtering
# ==================================================================================================================================================

# FILTERS:
# variant <50% missing data
# minQ 30
# individual with to much missing data  >0.5
# min DP 2
# bi-allelic
# variant <80% missing data
# mean DP between combined:3-13.568
# strand bias p-value combined< 0.01
# variant <85% missing data
# imiss >30%
# mac 5
# LD filter


# ----------------------------
# ------- concat.sh ----------
# ----------------------------

source $GDCstack
module load bcftools
module load vcftools/0.1.16-tc6l6nq

VCF_FILE="ac_ae_raw.vcf.gz"
# 51559526 SNPs

# make list of bcfs to combine
ls /SNPs/*.bcf > bcf_list

# Concatenate bcf files
bcftools concat --file-list bcf_list -Oz -o $VCF_FILE

echo "Finished combining. Start calculating stats $(date)"

# missingness per individual
vcftools --gzvcf ${VCF_FILE} --missing-indv --out ${VCF_FILE/.vcf.gz/} &

# mean coverage across the sites.
vcftools --gzvcf $VCF_FILE --site-mean-depth --out ${VCF_FILE/.vcf.gz/} &

# individual depth
vcftools --gzvcf ${VCF_FILE} --depth --out ${VCF_FILE/.vcf.gz/} &

# Wait for all background processes to complete
wait

echo "All vcftools analyses are done!"


# ----------------------------
# ------ interactive ---------
# ----------------------------

# plot missingnes and depth
# filter site missingness, mapping quality, minimum depth, biallelic sites and remove samples with missingness above cutoff
source $GDCstack
module load r


VCF_FILE="ac_ae_raw"

# plot individual missingness as a histogram
# Ind missing
awk 'NR > 1 {print $5}' "${VCF_FILE/.vcf.gz/}.imiss" | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_ind_miss.pdf")); 
    hist(data, main="Histogram of Missingness per individual", xlab="Missingness", col="steelblue", border="black", breaks=100); 
    dev.off()' "${VCF_FILE/.vcf.gz/}" &

# calculate average missingness
awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' ${VCF_FILE}.imiss
# 0.0893765
# list samples which have more than 50% missingness
tail +2 ${VCF_FILE}.imiss | awk '$5>0.5' | cut -f 1 > samples_imiss5
# Number of samples to be removed
wc -l samples_imiss5
# 3 


# DEPTH

# individual depth
awk 'NR > 1 {sum += $3; count++; 
  if (count == 1) {
    max = $3;
    min = $3;
  }
  if ($3 > max) max = $3; 
  if ($3 < min) min = $3; 
} END { 
  if (count > 0) {
    print "Mean Depth:", sum / count;
    print "Max Depth:", max;
    print "Min Depth:", min;
  }
}' ${VCF_FILE}.idepth
# Mean Depth: 4.59095
# Max Depth: 9.39831
# Min Depth: 0.0251447

tail -n+2 ${VCF_FILE}.idepth | awk -F '\t' '$3<5' | wc -l
# 09 individuals have an average depth below 5


# Ind depth
awk 'NR > 1 {print $3}' ${VCF_FILE/.vcf.gz/}.idepth | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_ind_dp.pdf")); 
    hist(data, main="Histogram of Mean Depth per individual", xlab="Mean Depth", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &


# mean coverage across sites
cut -f3 ${VCF_FILE/.vcf.gz/}.ldepth.mean | tail +2 > rawdepthpersite
# calculate quantiles
Rscript -e 'quantile (as.numeric (readLines ("stdin")), probs = c(0.01, 0.05, 0.5, 0.95, 0.99))' < rawdepthpersite
#        1%        5%       50%       95%       99% 
# 0.055336  0.335968  3.766800  7.122530 12.221300 

awk 'NR > 1 {print $3}' ${VCF_FILE/.vcf.gz/}.ldepth.mean | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_site_dp.pdf"));
    hist(data, main="Histogram of Mean Depth per Site", xlab="Mean Depth", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &

# ---------------------------

# ----------------------------
# ------- filt_1.sh ----------
# ----------------------------

# filter for site missingness, mapping quality, minimum depth, biallelic sites and remove samples with missingness above cutoff

VCF_FILE="ac_ae_raw.vcf.gz"

# filter for site missingness, quality, minimum depth, biallelic sites and remove samples with missingness above cutoff
vcftools --gzvcf $VCF_FILE --max-missing 0.5 --minQ 30 --minDP 2 --max-alleles 2 --min-alleles 2 --remove samples_imiss5 --recode --recode-INFO-all --stdout | bgzip -c > ${VCF_FILE/.vcf.gz/_miss5_Q30_DP2_bi_imiss5.vcf.gz}
# After filtering, kept 38349898 out of a possible 51559526 Sites 

# stats
VCF_FILE2=${VCF_FILE/.vcf.gz/_miss5_Q30_DP2_bi_imiss5.vcf.gz}

# mean coverage across the sites.
vcftools --gzvcf ${VCF_FILE2} --site-mean-depth --out ${VCF_FILE2/.vcf.gz/} &

# missingnes per site
vcftools --gzvcf ${VCF_FILE2} --missing-site --out ${VCF_FILE2/.vcf.gz/} &

wait


# 38349899 SNps left 


# ---------------------------

# ----------------------------
# ------ check output ---------
# ----------------------------


module load r

VCF_FILE="ac_ae_raw_miss5_Q30_DP2_bi_imiss5.vcf.gz"

## DEPTH ##

# mean coverage across sites
cut -f3 ${VCF_FILE/.vcf.gz/}.ldepth.mean | tail +2 > meandepthpersite

# make a histrogram (cuting x-axis of at 4*mean)
awk 'NR > 1 {print $3}' ${VCF_FILE/.vcf.gz/}.ldepth.mean | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_site_dp.pdf"));
    hist(data, main="Histogram of Mean Depth per Site", xlab="Mean Depth", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &


# calculate quantiles
Rscript -e 'quantile (as.numeric (readLines ("stdin")), probs = c(0.01, 0.05, 0.5, 0.95, 0.99))' < meandepthpersite
#     1%      5%     50%     95%     99% 
# 2.20000 3.09697 4.92727 6.70909 8.86667 

# calculate mean
Rscript -e 'mean (as.numeric (readLines ("stdin")))' < meandepthpersite
# 4.968937 - 3*mean = 14.9

#potential cuts off:
# 1% and 99% quantile
# min dp 2-3
# max dp 2-3* mean

# min dp=3 and max dp as 99% quantile
awk '$1 >= 3 && $1 <= 8.86667' meandepthpersite | wc -l
# 36285947



#Let's check the missingness across the sites.
tail +2 ${VCF_FILE/.vcf.gz/}.lmiss | cut -f6 > totalmissingsite
# calculate 0.5 and 0.95 quantiles


# plot histogram
# SNP  missing
awk 'NR > 1 {print $6}' ${VCF_FILE/.vcf.gz/}.lmiss | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_site_miss.pdf")); 
    hist(data, main="Histogram of Missingness per Site", xlab="Missingness", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &

Rscript -e 'quantile (as.numeric (readLines ("stdin")),probs = c(0.05, 0.1, 0.5, 0.90, 0.95))' < totalmissingsite
#       5%       10%       50%       90%       95% 
# 0.0424242 0.0545455 0.1030300 0.2060610 0.2848480 

awk '$1 <= 0.2' totalmissingsite | wc -l
awk '$1 <= 0.1' totalmissingsite | wc -l
# 34377322
# 18853524


# ---------------------------

# ----------------------------
# ------- filt_2.sh ----------
# ----------------------------

# filter site based on missingess and mean site depth between 0.05 and 0.95 quantiles
vcftools --gzvcf ${VCF_FILE} --max-missing 0.8 --min-meanDP 3 --max-meanDP 8.86667 --mac 2 --recode-INFO-all --recode --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss8_meanDP_mac2.vcf.gz} 
# After filtering, kept 33119745 out of a possible 38349898 Sites

# filter for stranbias
bcftools filter -e 'INFO/PV4[0]< 0.01' ${VCF_FILE/.vcf.gz/_miss8_meanDP_mac2.vcf.gz} -Oz -o ${VCF_FILE/.vcf.gz/_miss8_meanDP_mac2_SB.vcf.gz}

# ----------------------------
# -------- stats.sh ----------
# ----------------------------

module load vcftools

VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB.vcf.gz"


echo "Finished filtering and calculating individual stats $(date)"
vcftools --gzvcf $VCF_FILE --depth --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --missing-indv --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --het --out $OUT_DIR/$OUT_PREFIX &
#vcftools --gzvcf $VCF_FILE --relatedness --out $OUT_DIR/$OUT_PREFIX &
#vcftools --gzvcf $VCF_FILE --relatedness2 --out $OUT_DIR/$OUT_PREFIX &

echo "Calculating allele stats $(date)"
vcftools --gzvcf $VCF_FILE --counts2 --out $OUT_DIR/$OUT_PREFIX &

echo "Calculating site stats $(date)"
#vcftools --gzvcf $VCF_FILE --site-depth --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --site-mean-depth --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --missing-site --out $OUT_DIR/$OUT_PREFIX &
#vcftools --gzvcf $VCF_FILE --singletons --out $OUT_DIR/$OUT_PREFIX &
vcftools --gzvcf $VCF_FILE --hardy --out $OUT_DIR/$OUT_PREFIX &


tail -n+2 ${VCF_FILE/.vcf.gz/.imiss} | awk -F '\t' '$5>0.3' | wc -l

tail -n+2 ${VCF_FILE/.vcf.gz/.lmiss}  | awk -F '\t' '$6>0.05' | wc -l

tail -n+2 ${VCF_FILE/.vcf.gz/.frq} | awk -F '\t' '$6<0.01' | wc -l

tail -n+2 ${VCF_FILE/.vcf.gz/.frq.count} | awk -F '\t' '$6<5' | wc -l

# ----------------------------
# ------- Plot stats ----------
# ----------------------------


source $GDCstack
module load r

VCF_FILE="ac_ae_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB.vcf.gz"
#VCF_FILE="ac_ae_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB_miss85_imiss3_mac5.vcf.gz"
VCF_FILE="cygnea_mac5.vcf.gz"
VCF_FILE="exulcerata_mac2.vcf.gz"
awk 'NR > 1 {print $3}' ${VCF_FILE/.vcf.gz/}.ldepth.mean | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_site_dp.pdf"));
    hist(data, main="Histogram of Mean Depth per Site", xlab="Mean Depth", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &

# SNP  missing
awk 'NR > 1 {print $6}' ${VCF_FILE/.vcf.gz/}.lmiss | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_site_miss.pdf")); 
    hist(data, main="Histogram of Missingness per Site", xlab="Missingness", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &

# SNP heterozygosity
awk 'NR > 1 {split($3, a, "/"); sum = a[1] + a[2] + a[3]; print a[2] / sum}' ${VCF_FILE/.vcf.gz/}.hwe | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_site_het.pdf")); 
    hist(data, main="Histogram of Heterozygosity per Site", xlab="Heterozygosity", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &

# Allele count
awk 'NR > 1 {print $6}' ${VCF_FILE/.vcf.gz/}.frq.count | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_allele_count.pdf")); 
    hist(data, main="Histogram of Allele Counts", xlab="Counts", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &

# Allele frequency
awk 'NR > 1 {print $6/$4}' "${VCF_FILE/.vcf.gz/}.frq.count" | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_allele_frq.pdf")); 
    hist(data, main="Histogram of Allele Frequency", xlab="Frequency", col="steelblue", border="black", breaks=100); 
    dev.off()' "${VCF_FILE/.vcf.gz/}" &

# Ind missing

awk 'NR > 1 {print $5}' "${VCF_FILE/.vcf.gz/}.imiss" | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_ind_miss.pdf")); 
    hist(data, main="Histogram of Missingness per individual", xlab="Missingness", col="steelblue", border="black", breaks=100); 
    dev.off()' "${VCF_FILE/.vcf.gz/}" &

# Ind heterozygosity
awk 'NR > 1 {print 1 - ($2 / $4)}' ${VCF_FILE/.vcf.gz/}.het | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_ind_het.pdf")); 
    hist(data, main="Histogram of Heterozygosity per individual", xlab="Heterozygosity", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &

# Ind depth
awk 'NR > 1 {print $3}' ${VCF_FILE/.vcf.gz/}.idepth | \
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); 
    data <- as.numeric(readLines("stdin")); 
    pdf(paste0(args[1], "_ind_dp.pdf")); 
    hist(data, main="Histogram of Mean Depth per individual", xlab="Mean Depth", col="steelblue", border="black", breaks=100); 
    dev.off()' ${VCF_FILE/.vcf.gz/} &


# ----------------------------
# ------- filt_3.sh ----------

VCF_FILE="ac_ae_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB.vcf.gz"

# filter for missingness
vcftools --gzvcf ${VCF_FILE} --max-missing 0.85  --recode-INFO-all --recode --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss85.vcf.gz} 

echo "Finished filtering and calculating individual stats $(date)"
vcftools --gzvcf ${VCF_FILE/.vcf.gz/_miss85.vcf.gz} --missing-indv --out ${VCF_FILE/.vcf.gz/_miss85} 

tail +2 ${VCF_FILE/.vcf.gz/_miss85}.imiss | awk '$5>0.3' | cut -f 1 > samples_imiss3

vcftools --gzvcf ${VCF_FILE/.vcf.gz/_miss85.vcf.gz} --remove samples_imiss3 --recode --recode-INFO-all --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss85_imiss3_mac5.vcf.gz} 



# ------------------------------
# ------- subsampling ----------
# ----------------------------

# make 
# - bed files
# - smaller vcfs 
# - species subest

source $GDCstack
module load plink
module load vcftools/0.1.16-tc6l6nq

VCF_FILE="ac_ae_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB_miss85_imiss3_mac5.vcf.gz"
BED_FILE="cygnea_exulcerata"

#  make a bed file
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out ${OUT_DIR}/${BED_FILE}

# make a vcf file which has less info field and thus is smaller and easier to read into memory
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --recode vcf-iid --out ${OUT_DIR}/$BED_FILE

###### CYGNEA #########

awk -F',' '$1 == "Anodonta cygnea" {print $NF "\t" $NF}' metadata_all.csv > cygnea_samples
plink --bfile ${OUT_DIR}/$BED_FILE --allow-extra-chr --keep cygnea_samples --make-bed --pca --out ${OUT_DIR}/cygnea

grep -vE "AC83|AC90" cygnea_samples > cygnea_samples_exc_hy
plink --bfile ${OUT_DIR}/$BED_FILE --allow-extra-chr --keep cygnea_samples_exc_hy --mac 5 --make-bed --pca --recode vcf-iid --out ${OUT_DIR}/cygnea_mac5 


###### EXULCERATA #########

awk -F',' '$1 == "Anodonta exulcerata" {print $NF "\t" $NF}' metadata_all.csv  > exulcerata_samples
plink --bfile ${OUT_DIR}/$BED_FILE --allow-extra-chr --keep exulcerata_samples --make-bed --pca --out ${OUT_DIR}/exulcerata 
plink --bfile ${OUT_DIR}/$BED_FILE --allow-extra-chr --keep exulcerata_samples --mac 2 --make-bed --pca --recode vcf-iid  --out ${OUT_DIR}/exulcerata_mac2 

###### EXULCERATA & CYGNEA #########


# ----------------------------
# -------- pcaone_ld.sh ----------
# ----------------------------

# filter for LD but account for population (species) structure

source $GDCstack
module load plink


# split vcf file in two for memory usage
BED_FILE="cygnea_exulcerata"


$HOME/software/PCAone/PCAone -b $BED_FILE -k 3 -n 1 --ld-r2 0.1 --ld-bp 100000 --out $BED_FILE

awk '{print $2 }' ${BED_FILE}.ld.prune.out > snp2rm

# Finally filter out SNPs in LD
#  make a bed file with ld filter ( here I also do pca )
plink --bfile $BED_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --exclude snp2rm --make-bed --pca --out ${BED_FILE}_ld

# make vcf file from filtered bed file
plink --bfile --out ${BED_FILE}_ld --allow-extra-chr --recode vcf-iid --out ${BED_FILE}_ld

BED_FILE="cygnea_mac5"
# make vcf file from filtered bed file
plink --bfile $BED_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --exclude snp2rm --make-bed --pca --out ${BED_FILE}_ld

BED_FILE="exulcerata_mac2"
# make vcf file from filtered bed file
plink --bfile $BED_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --exclude snp2rm --make-bed --pca --out ${BED_FILE}_ld









# ==================================================================================================================================================
# 7. Analysis
# ==================================================================================================================================================



# 0. Quality (LD decay, allele frequency)
# 1. Population structure (PCA, admixture, Fst, RAxML)
# 2. Genetic diversity (dxy, pi, Ho, He, Polymorphic sites)
# 3. Inbreeding (FIS, ROH, genetic load)
# 4. Effective population size (currentNe)



# Eg of how to chain jobs
JOBID1=$(sbatch admix.sh cygnea_exulcerata_ld | awk '{print $4}')
JOBID2=$(sbatch --dependency=afterany:$JOBID1 admix.sh cygnea_mac5_ld | awk '{print $4}')
JOBID3=$(sbatch --dependency=afterany:$JOBID2 admix.sh | awk '{print $4}')
JOBID4=$(sbatch --dependency=afterany:$JOBID3 admix.sh | awk '{print $4}')
# Other Dependency Types
# afterok: Job starts if the specified job completes successfully.
# afterany: Job starts after the specified job finishes, regardless of exit status.
# afternotok: Job starts if the specified job fails.
# singleton: Ensures only one instance of a job with the same name runs at a time.

# ==================================================================================================================================================
# Quality
# ==================================================================================================================================================

# ----------------------------
# ------------AF-------------- DONE
# ----------------------------
out="afreq"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink2

BED_FILE="cygnea_exulcerata"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out ${out}/$BED_FILE

BED_FILE="cygnea_mac5"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out ${out}/$BED_FILE

BED_FILE="exulcerata_mac2"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out ${out}/$BED_FILE

# ----------------------------
# --------LD decay------------ DONE
# ----------------------------

out="ld"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink

BED_FILE="cygnea_exulcerata"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.25 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --make-bed --out ${out}/${BED_FILE}_maf_thin
awk '{print $7, $5-$2}' ${out}/${BED_FILE}_maf_thin.ld > ${BED_FILE}_maf_thin_ld_with_distance.txt 

BED_FILE="cygnea_mac5"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.25 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --make-bed --out ${out}/${BED_FILE}_maf_thin
awk '{print $7, $5-$2}' ${out}/${BED_FILE}_maf_thin.ld > ${BED_FILE}_maf_thin_ld_with_distance.txt 

BED_FILE="exulcerata_mac2"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.25 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --make-bed --out ${out}/${BED_FILE}_maf_thin
awk '{print $7, $5-$2}' ${out}/${BED_FILE}_maf_thin.ld > ${out}/${BED_FILE}_maf_thin_ld_with_distance.txt 

# This command creates a files with
# r² value (R2)
# distance between POS_B and POS_A

# in case the output is still to large to plot or process you can subset like this

awk 'NR==1{print; next} rand() <= 500000000/NR' ${BED_FILE}_maf_thin_ld_with_distance.txt  > ${BED_FILE}_maf_thin_ld_with_distance_subset.txt 


# ----------------------------
# ----Allele DP + Het frequency----- DONE
# ----------------------------

source $GDCstack
module load bcftools

VCF_FILE=ac_ae_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB_miss85_imiss3_mac5
tabix -p vcf ${VCF_FILE}.vcf.gz
# subset vcf file
bcftools view -H ${VCF_FILE}.vcf.gz | shuf -n 50000000 | cut -f1,2 > ${VCF_FILE}_5Msnsp.txt &


# can be necessary if the metadala_all was made in excel
dos2uniz metadata_all.csv

awk -F',' '$1 == "Anodonta cygnea" && $9 == "passed" {print $NF}' metadata_all.csv > cygnea_samples

grep -vE "AC83|AC90" cygnea_samples > cygnea_samples_exc_hy

awk -F',' '$1 == "Anodonta exulcerata"&& $9 == "passed" {print $NF}' metadata_all.csv  > exulcerata_samples



bcftools query -S exulcerata_samples -T ${VCF_FILE}_5Msnsp.txt  -f '%CHROM:%POS\t[%GT:%AD;]\n' ${VCF_FILE}.vcf.gz > ae_allele_depth.txt
bcftools query -S cygnea_samples_exc_hy -T ${VCF_FILE}_5Msnsp.txt  -f '%CHROM:%POS\t[%GT:%AD;]\n' ${VCF_FILE}.vcf.gz > ac_allele_depth.txt



awk '
BEGIN {
    print "SNP\tn_01\tDP_0\tDP_1\tn_total"
}
{
    split($2, samples, ";");
    sum_01 = 0;
    count_01 = 0;
    sum_01_allele1 = 0;
    sum_01_allele0 = 0;
    n_total = 0;
    
    for (i in samples) {
        split(samples[i], info, ":");
        genotype = info[1];
        allele_depth = info[2];
        split(allele_depth, depths, ",");
        depth_0 = depths[1];  # Depth for allele 0
        depth_1 = depths[2];  # Depth for allele 1

        # Count all genotyped individuals (exclude missing "./.")
        if (genotype != "./.") {
            n_total++;
        }

        # Count heterozygotes (0/1)
        if (genotype == "0/1") {
            sum_01_allele0 += depth_0;
            sum_01_allele1 += depth_1;
            count_01++;
        }
    }

    # Print results for each SNP
    print $1 "\t" count_01 "\t" sum_01_allele0 "\t" sum_01_allele1 "\t" n_total;
}' ac_allele_depth.txt > ac_ad.txt



awk '
BEGIN {
    print "SNP\tn_01\tDP_0\tDP_1\tn_total"
}
{
    split($2, samples, ";");
    sum_01 = 0;
    count_01 = 0;
    sum_01_allele1 = 0;
    sum_01_allele0 = 0;
    n_total = 0;
    
    for (i in samples) {
        split(samples[i], info, ":");
        genotype = info[1];
        allele_depth = info[2];
        split(allele_depth, depths, ",");
        depth_0 = depths[1];  # Depth for allele 0
        depth_1 = depths[2];  # Depth for allele 1

        # Count all genotyped individuals (exclude missing "./.")
        if (genotype != "./.") {
            n_total++;
        }

        # Count heterozygotes (0/1)
        if (genotype == "0/1") {
            sum_01_allele0 += depth_0;
            sum_01_allele1 += depth_1;
            count_01++;
        }
    }

    # Print results for each SNP
    print $1 "\t" count_01 "\t" sum_01_allele0 "\t" sum_01_allele1 "\t" n_total;
}' ae_allele_depth.txt > ae_ad.txt





# heterozygosity per genotype aa ab bb
# allele depth difference in heterozygots
# anatina - raxml
# cygnea exulcerata -r raxml


# ==================================================================================================================================================
# 1 Population structure 
# ==================================================================================================================================================


# ----------------------------
# ------------PCA------------- DONE
# ----------------------------
out="pca"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink

BED_FILE="cygnea_exulcerata_ld"
BED_FILE="cygnea_mac5_ld"
BED_FILE="exulcerata_mac2_ld"

# create pca and files that allow proper percentage variance calculations
plink --bfile $BED_FILE --allow-extra-chr --pca var-wts --make-rel --out ${out}/$BED_FILE
#plink --bfile $BED_FILE --allow-extra-chr --pca --out $BED_FILE

#Calculate sum of variance from relation covariance matrix:
sum1=$(awk '{sum+=$NR;}END{print sum}' ${out}/${BED_FILE}.rel)

#Calculate percentage variance explained (pve) and write to file
while read line; do
    echo "$(echo "scale=4; $line / $sum1 * 100" | bc)"
done < ${out}/${BED_FILE}.eigenval > ${out}/${BED_FILE}.pve

# calculate missingness to be plotted on PCs
plink --bfile $BED_FILE --allow-extra-chr --missing --out ${out}/$BED_FILE 

# ----------------------------
# ---------ADMIXTURE---------- DONE
# ----------------------------

source $GDCstack
module load plink

BED_FILE="cygnea_exulcerata_ld"
BED_FILE="cygnea_mac5_ld"

# make folder and move files

out="admix"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi


cp -v ${BED_FILE}.* ./admix/
cd admix

# ADMIXTURE does not accept chromosome names that are not human chromosomes. 
# We will thus just exchange the first column by 0

awk '{$1="0";print $0}' ${BED_FILE}.bim > ${BED_FILE}.bim.tmp
mv ${BED_FILE}.bim.tmp ${BED_FILE}.bim


# test run
#source $GDCstack
#module load admixture
#admixture --cv=5 -j8 ${BED_FILE}.bed 2 

BED_FILE="cygnea_exulcerata_ld"
JOBID1=$(sbatch admix.sh $BED_FILE | awk '{print $4}')

BED_FILE="cygnea_mac5_ld"
JOBID2=$(sbatch --dependency=afterany:$JOBID1 admix.sh $BED_FILE | awk '{print $4}')

# note that admixture needs a lot of memory for high Ks, especially when calculating cv for many SNPs
# if it keeps running out of memory try thinning the data, either by LD or by position.

# get cross validation values

awk ' /CV/ {print $3,$4}' *${BED_FILE}*out | sed -e 's/(K=//;s/)://' | sort -n -k 1 > ${BED_FILE}.cv


# ----------------------------
# ---------    FST  ---------- DONE
# ----------------------------

out="fst"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink2


BED_FILE="cygnea_exulcerata_ld"
BED_FILE="cygnea_exulcerata"


# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 && $9 == "passed" {print $NF "\t" $NF "\t" $2 "_" $5}' metadata_all.csv > pheno

plink2 --bfile ${BED_FILE} --double-id --allow-extra-chr --pheno pheno --fst Pop --out ${out}/${BED_FILE}_pop &

# only do pairwise fst for pops with 5 or more individuals individual
awk -F, '
BEGIN { OFS="\t"; print "FID", "IID", "Pop" } 
NR > 1 && $9 == "passed" { count[$2 "_" $5]++; data[NR] = $NF "\t" $NF "\t" $2 "_" $5 }
END {
  for (i in data) {
    split(data[i], fields, "\t")
    if (count[fields[3]] >= 5)
      print data[i]
  }
}' metadata_all.csv > pheno_5
plink2 --bfile ${BED_FILE} --double-id --allow-extra-chr --pheno pheno_5 --fst Pop --out ${out}/${BED_FILE}_pop5 


# ----------------------------
# ---------RAxML---------- DONE
# ----------------------------

# https://github.com/amkozlov/raxml-ng/wiki/Tutorial 

out="tree"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi


BED_FILE="cygnea_exulcerata_ld"


source $GDCstack
module load raxml-ng/1.2.2 
module load python
module load plink
module load bcftools

plink --bfile $BED_FILE --allow-extra-chr --recode vcf-iid --out ${BED_FILE}

#to leave only SNPs that were present at least once in our dataset as  homozygous for the reference allele, and homozygous for the alternative allele, as required by RAxML.
bcftools view -i 'COUNT(GT="RR")>0 & COUNT(GT="AA")>0' ${BED_FILE}.vcf -Ov -o ${BED_FILE}_filtered.vcf  # Uncompressed VCF
bcftools view -H ${BED_FILE}_filtered.vcf | wc -l 
# 766190


#cat $BED_FILE.bim | cut -f 2 > sites
#shuf sites -n 100000 | sort > sites100000
#plink --bfile $BED_FILE --allow-extra-chr --extract sites100000 --recode vcf-iid --out ${BED_FILE}_100k


BED_FILE="cygnea_exulcerata_ld_filtered"
PHY_FILE="cygnea_exulcerata_ld_filtered.min4.phy"

# convert vcf format to phylip for our raxml analyses.
python $HOME/software/vcf2phylip-master/vcf2phylip.py -i ${BED_FILE}.vcf --output-folder ${out}


# check that the MSA can actually be read
raxml-ng --check --msa $PHY_FILE --model GTR+ASC_LEWIS --prefix T1

# Compress alignment patterns and store MSA in the binary format (RAxML Binary Alignment, RBA
# Estimate memory requirements and optimal number of CPUs/threads
raxml-ng --parse --msa $PHY_FILE --model GTR+ASC_LEWIS --prefix T2

raxml-ng --msa T2.raxml.rba --model GTR+ASC_LEWIS --prefix T3 --threads 22 --seed 2 # takes about 3mi


# ----------------------------
# --------- IBS/Hamming  ----------
# ----------------------------


BED_FILE="cygnea_exulcerata_ld"


source $GDCstack
module load plink


plink --bfile $BED_FILE --allow-extra-chr   --distance square --out ${BED_FILE}

# ==================================================================================================================================================
# 2. Diversity
# ==================================================================================================================================================


# ----------------------------
# ------------dxy------------- SKIP
# ----------------------------



# --------------------------------
# ---nucleotide diversity -------- DONE
# --------------------------------
# pi.sh


out="pi"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load vcftools/0.1.16-tc6l6nq

BED_FILE="cygnea_exulcerata"

# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 && $9 == "passed" {print $NF "\t" $NF "\t" $2 "_" $5}' metadata_all.csv > pheno
# Extract unique groups (populations) from the third column of the pheno file
tail -n +2 pheno | cut -f 3 | sort | uniq > pops

grep -v "AA_" pops > pops_AC_AE


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
    fi
done < pops_AC_AE


# concatenate the files
echo -e "pop\tCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI" > ${out}/${BED_FILE}_pi_10kb.txt
while read -r group; do
tail -n +2 ${out}/pi_10kb_${group}.windowed.pi | awk -v grp="$group" '{print grp "\t" $0}' >> ${out}/${BED_FILE}_pi_10kb.txt
done < pops_AC_AE


# ----------------------------
# -----------Het-------------- DONE
# ----------------------------

# fis.sh
# same as below FIS




# ==================================================================================================================================================
# 3. Inbreeding
# ==================================================================================================================================================


# ----------------------------
# -----------FIS-------------- DONE
# ----------------------------
# fis.sh

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


# ----------------------------
# ----------- ROH ------------ DONE
# ----------------------------

out="roh"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load bcftools
module load vcftools/0.1.16-tc6l6nq


VCF_FILE="cygnea_mac5.vcf"
VCF_FILE="exulcerata_mac2.vcf"

bgzip $VCF_FILE
tabix -p vcf ${VCF_FILE}.gz

# add allele frequency to file
bcftools +fill-tags ${VCF_FILE}.gz -Oz -o ${VCF_FILE/.vcf/}_AF.vcf.gz -- -t all




# This limits ROH calculations to sites with 30 or more genotyped individuals.
# use the AF field for allele frequency
# only iutput RG
bcftools roh --GTs-only 30 --AF-tag AF ${VCF_FILE/.vcf/}_AF.vcf.gz --output-type r -o ${out}/${VCF_FILE}_roh.txt


# ----------------------------
# -------relatedness --------- DONE
# ----------------------------
out="kinship"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink2
BED_FILE="cygnea_mac5_ld"
BED_FILE="exulcerata_mac2_ld"
BED_FILE="cygnea_exulcerata_ld"
BED_FILE="cygnea_exulcerata"
BED_FILE="cygnea_mac5"
BED_FILE="exulcerata_mac2"


# KING kinship estimator
plink2 --bfile $BED_FILE --allow-extra-chr --make-king-table counts --make-king square --out ${out}/$BED_FILE

# to check  for high values
for file in ${out}/*.kin0; do
    count=$(awk '$NF > 0.2' "$file" | wc -l)
    echo "$file: $count"
done



# ==================================================================================================================================================
# 4. Effective population size
# ==================================================================================================================================================

# ----------------------------
# -------- currentNe --------- DONE
# ----------------------------


out="ne"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

# subset vcf file
bcftools view -H ${VCF_FILE}.vcf | shuf -n 1999999 | cut -f1,2 > ${out}/${VCF_FILE}_2Msnsp.txt

# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 && $9 == "passed" {print $NF "\t" $NF "\t" $2 "_" $5}' metadata_all.csv > pheno

# Extract unique groups (populations) from the third column of the pheno file
tail -n +2 pheno | cut -f 3 | sort | uniq -c | awk '$1 >= 10 {print $2}' > pops_10

grep "AC_" pops_10 > pops_10_AC
grep "AE_" pops_10 > pops_10_AE

VCF_FILE="cygnea_mac5"
VCF_FILE="exulcerata_mac2"


# ne.sh

VCF_FILE="cygnea_mac5"
pops=pops_10_AC
out=ne


# how to reference each sample by array number
IDX=$SLURM_ARRAY_TASK_ID
group=`sed -n ${IDX}p < $pops`

echo "Processing file: $pops"

source $GDCstack
module load bcftools

# Create the 'keep' file for the group
awk -v grp="$group" '$3 == grp {print $1}' pheno > ${out}/keep$group
    
bcftools view -T ${out}/${VCF_FILE}_2Msnsp.txt -S ${out}/keep${group} -o ${out}/${group}.vcf -O v ${VCF_FILE}.vcf.gz

echo "Staring Ne calculations of $pops"

$HOME/software/currentNe/currentNe -t 8 ${out}/${group}.vcf 19



# Once finished runs this on output

# To see Ne estimate an CI
grep -A9 "# Ne point estimate:" *_currentNe_OUTPUT.txt

# To only see Ne estimate
grep -A1 "# Ne point estimate:" *_currentNe_OUTPUT.txt

# Ne estimates without genetic map
for group in $(< pops_10_AC); do
    ne=$(awk 'NR==50' "${group}_currentNe_OUTPUT.txt")
    CI_10=$(awk 'NR==56' "${group}_currentNe_OUTPUT.txt")
    CI_90=$(awk 'NR==58' "${group}_currentNe_OUTPUT.txt")
    echo "$group $ne CI:${CI_10}-${CI_90}"
done

# AC_ETD 6.24 CI:5.01-7.89
# AC_ETR 4.00 CI:5.01-5.01
# AC_ROT 21.52 CI:16.61-27.88
# AC_SAM 4.00 CI:5.01-5.01
# AC_UFS 23.66 CI:19.52-28.69
# AC_VIL 4.00 CI:5.01-5.01
# AE_AGN 7.41 CI:6.43-8.54


# ==================================================================================================================================================
# TEST ZONE
# ==================================================================================================================================================

source $GDCstack
module load bcftools

VCF_FILE=ac_ae_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB_miss85_imiss3_mac5
tabix -p vcf ${VCF_FILE}.vcf.gz
# subset vcf file
bcftools view -H ${VCF_FILE}.vcf.gz | shuf -n 50000000 | cut -f1,2 > ${VCF_FILE}_5Msnsp.txt &


# can be necessary if the metadala_all was made in excel
dos2uniz metadata_all.csv

awk -F',' '$1 == "Anodonta cygnea" && $9 == "passed" {print $NF}' metadata_all.csv > cygnea_samples

grep -vE "AC83|AC90" cygnea_samples > cygnea_samples_exc_hy

awk -F',' '$1 == "Anodonta exulcerata"&& $9 == "passed" {print $NF}' metadata_all.csv  > exulcerata_samples



bcftools query -S exulcerata_samples -T ${VCF_FILE}_5Msnsp.txt  -f '%CHROM:%POS\t[%GT:%AD;]\n' ${VCF_FILE}.vcf.gz > ae_allele_depth.txt
bcftools query -S cygnea_samples_exc_hy -T ${VCF_FILE}_5Msnsp.txt  -f '%CHROM:%POS\t[%GT:%AD;]\n' ${VCF_FILE}.vcf.gz > ac_allele_depth.txt



awk '
BEGIN {
    print "SNP\tn_01\tDP_0\tDP_1\tn_total"
}
{
    split($2, samples, ";");
    sum_01 = 0;
    count_01 = 0;
    sum_01_allele1 = 0;
    sum_01_allele0 = 0;
    n_total = 0;
    
    for (i in samples) {
        split(samples[i], info, ":");
        genotype = info[1];
        allele_depth = info[2];
        split(allele_depth, depths, ",");
        depth_0 = depths[1];  # Depth for allele 0
        depth_1 = depths[2];  # Depth for allele 1

        # Count all genotyped individuals (exclude missing "./.")
        if (genotype != "./.") {
            n_total++;
        }

        # Count heterozygotes (0/1)
        if (genotype == "0/1") {
            sum_01_allele0 += depth_0;
            sum_01_allele1 += depth_1;
            count_01++;
        }
    }

    # Print results for each SNP
    print $1 "\t" count_01 "\t" sum_01_allele0 "\t" sum_01_allele1 "\t" n_total;
}' ac_allele_depth.txt > ac_ad.txt



awk '
BEGIN {
    print "SNP\tn_01\tDP_0\tDP_1\tn_total"
}
{
    split($2, samples, ";");
    sum_01 = 0;
    count_01 = 0;
    sum_01_allele1 = 0;
    sum_01_allele0 = 0;
    n_total = 0;
    
    for (i in samples) {
        split(samples[i], info, ":");
        genotype = info[1];
        allele_depth = info[2];
        split(allele_depth, depths, ",");
        depth_0 = depths[1];  # Depth for allele 0
        depth_1 = depths[2];  # Depth for allele 1

        # Count all genotyped individuals (exclude missing "./.")
        if (genotype != "./.") {
            n_total++;
        }

        # Count heterozygotes (0/1)
        if (genotype == "0/1") {
            sum_01_allele0 += depth_0;
            sum_01_allele1 += depth_1;
            count_01++;
        }
    }

    # Print results for each SNP
    print $1 "\t" count_01 "\t" sum_01_allele0 "\t" sum_01_allele1 "\t" n_total;
}' ae_allele_depth.txt > ae_ad.txt





# heterozygosity per genotype aa ab bb
# allele depth difference in heterozygots
# anatina - raxml
# cygnea exulcerata -r raxml






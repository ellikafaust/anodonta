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
source /cluster/project/gdc/shared/stack/GDCstack.sh
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

module load stack
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
#SBATCH --job-name=index 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G
#SBATCH --time=1:00:00
#SBATCH --output=index.o%A
#SBATCH --error=index.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julie.conrads@eawag.ch

module load bwa-mem2
bwa-mem2 index /cluster/scratch/jconrads/Aanat_demux_hifi_l3.bp.p_ctg.fa
'> index.sh

sbatch index.sh
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
ls *_R1_001.fastq.gz |sed "s/_R1_001.fastq.gz//" > sample_list

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
ls *_R1_001.fastq.gz |sed "s/_R1_001.fastq.gz//" > sample_list

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
module load multiqc

multiqc mapping/stats/ -n multiqc
multiqc mapping/stats_dup/ -n multiqc_dup
multiqc mapping/statsQ20/ -n multiqc_nodup_Q20


# --------  tar archive and copy mapping stats, logs and multiqc files to project ----------



# ==================================================================================================================================================
# 5. Variant calling
# ==================================================================================================================================================



# --------  copy bam files to scratch (will take some time) ----------

index=/cluster/work/gdc/shared/p976/index/Aanat_demux_hifi_l3.bp.p_ctg.fa.fai

# prepare sorted scaffold list from index
cut -f 1,2 *.fai | sort -k2 -n -r > sorted_scaffolds
# prepare list of bam files
ls bams/*bam > bam_list

# if there are more bams than the once I am interested in
awk -F',' '$1 == "Anodonta anatina" {print $NF}' groups_aa.csv | grep -Ff - bam_list > anatina_bams

# split sorted scaffolds according to size
head -n 100 sorted_scaffolds | cut -f 1 > largest_scaffolds_100
tail -n +101 sorted_scaffolds | cut -f 1 | shuf  > smaller_scaffolds
split -l 100 smaller_scaffolds scaffold_chunk_
# make the names better
counter=1
for file in scaffold_chunk_*; do
    mv "$file" "scaffold_chunk_${counter}"
    ((counter++))
done

# run bcfcalls.sh with largest 50 scaffolds
# run bcfcall_loop.sh with remaining scaffolds

# bcfcalls.sh will do the follwing
# 1. generate genotype likelihoods with bcftools mpileup
#		skipping indels, with a max read-depth of 30
# 2. call genotypes with bcftools call
#		using multiallelic-caller method, only outputting variant sites
sbatch bcfcalls.sh
sbatch bcfcalls_loop.sh

# after completed run check all log files for potential error messages and samples which did not complete
grep "Job: " *.out

# --------  copy raw bcf files to project ----------



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
# concatenate bcf files


VCF_FILE="aa_raw.vcf.gz"
# Concatenate bcf files
bcftools concat --file-list bcf_list -Oz -o $VCF_FILE

# missingness per individual
vcftools --gzvcf ${VCF_FILE} --missing-indv --out ${VCF_FILE/.vcf.gz/} &

# mean coverage across the sites.
vcftools --gzvcf $VCF_FILE --site-mean-depth --out ${VCF_FILE/.vcf.gz/} &

# individual depth
vcftools --gzvcf ${VCF_FILE} --depth --out ${VCF_FILE/.vcf.gz/} &

wait
# Optional

#parallel -j 3 vcftools --gzvcf ${VCF_FILE} {} --out ${VCF_FILE} ::: \
#    --missing-indv \
#    --site-mean-depth \
#    --depth
# ---------------------------

# ----------------------------
# ------ check output  ---------
# ----------------------------

# filter for site missingness, mapping quality, minimum depth, biallelic sites and remove samples with missingness above cutoff

VCF_FILE="aa_raw"

# plot individual missingness as a histogram
mawk '!/IN/' ${VCF_FILE}.imiss | cut -f5 > totalmissing

Rscript -e 'data <- as.numeric(readLines("stdin")); 
    pdf("aa_raw_miss_ind_hist.pdf"); 
    hist(data, main="Histogram of Missingness per individual", xlab="Frequency missingness", col="steelblue", border="black", breaks=100); 
    dev.off()' < totalmissing

# calculate average missingness
awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' ${VCF_FILE}.imiss
# 0.22

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
# Mean Depth: 3.70406
# Max Depth: 5.98197
# Min Depth: 0.28318

tail -n+2 ${VCF_FILE}.idepth | awk -F '\t' '$3<5' | wc -l
# 237 individuals have an average depth below 5

tail -n+2 ${VCF_FILE}.idepth | cut -f3 > inddepth
Rscript -e 'data <- as.numeric(readLines("stdin")); 
    pdf("aa_raw_depth_ind_hist.pdf"); 
    hist(data, main="Histogram of Mean Depth per individual", xlab="Mean Depth", col="steelblue", border="black", breaks=100); 
    dev.off()' < inddepth

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load r

# mean coverage across sites
cut -f3 ${VCF_FILE}.ldepth.mean | tail +2 > rawdepthpersite
# calculate quantiles
Rscript -e 'quantile (as.numeric (readLines ("stdin")), probs = c(0.01, 0.05, 0.5, 0.95, 0.99))' < rawdepthpersite
#        1%        5%       50%       95%       99% 
# 0.055336  0.335968  3.766800  7.122530 12.221300 

Rscript -e 'data <- as.numeric(readLines("stdin")); 
    mean_depth <- mean(data); pdf("aa_raw_meandp_site_hist.pdf"); 
    hist(data, main="Histogram of Mean Depth per Site", xlab="Mean Depth", col="steelblue", border="black", breaks=100, xlim=c(0, 4*mean_depth)); 
    dev.off()' < rawdepthpersite

# ---------------------------

# ----------------------------
# ------- filt_1.sh ----------
# ----------------------------

# filter for site missingness, mapping quality, minimum depth, biallelic sites and remove samples with missingness above cutoff

VCF_FILE="aa_raw.vcf.gz"

# filter for site missingness, mapping quality, minimum depth, biallelic sites and remove samples with missingness above cutoff
vcftools --gzvcf $VCF_FILE --max-missing 0.5 --minQ 30 --minDP 2 --max-alleles 2 --min-alleles 2 --remove samples_imiss5 --recode --recode-INFO-all --stdout | bgzip -c > ${VCF_FILE/.vcf.gz/_miss5_Q30_DP2_bi_imiss5.vcf.gz}

# stats
VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5.vcf.gz"
# calculate missingness per individual
vcftools --gzvcf ${VCF_FILE} --missing-site --out ${VCF_FILE/.vcf.gz/} &

# mean coverage across the sites.
vcftools --gzvcf ${VCF_FILE} --site-mean-depth --out ${VCF_FILE/.vcf.gz/} &

wait


# After filtering, kept 250 out of 253 Individuals
# Outputting VCF file...
# After filtering, kept 74290801 out of a possible 128768281 Sites

# ---------------------------

# ----------------------------
# ------ check output ---------

module load r

VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5.vcf.gz"

## DEPTH ##

# mean coverage across sites
cut -f3 ${VCF_FILE/.vcf.gz/}.ldepth.mean | tail +2 > meandepthpersite
# calculate quantiles
Rscript -e 'quantile (as.numeric (readLines ("stdin")), probs = c(0.01, 0.05, 0.5, 0.95, 0.99))' < meandepthpersite
#    1%     5%    50%    95%    99% 
# 1.960  2.312  4.728  7.908 13.568 

# calculate mean
Rscript -e 'mean (as.numeric (readLines ("stdin")))' < meandepthpersite
# anatina 4.852593 - 3*mean = 14.557779

# calculate sd
Rscript -e 'sd (as.numeric (readLines ("stdin")))' < meandepthpersite
# anatina 2.30699

#potential cuts off:
# 1% and 99% quantile
# min dp 2-3
# max dp 2-3* mean

# min dp=3 and max dp as 99% quantile
awk '$1 >= 3 && $1 <= 13.568' meandepthpersite | wc -l

# make a histrogram (cuting x-axis of at 4*mean)
Rscript -e 'data <- as.numeric(readLines("stdin")); 
    mean_depth <- mean(data); pdf("aa_raw_miss5_Q30_DP2_bi_imiss5_site_dp.pdf"); 
    hist(data, main="Histogram of Mean Depth per Site", xlab="Mean Depth", col="steelblue", border="black", breaks=100, xlim=c(0, 4*mean_depth)); 
    dev.off()' < meandepthpersite


#Let's check the missingness across the sites.
tail +2 ${VCF_FILE/.vcf.gz/}.lmiss | cut -f6 > totalmissingsite
# calculate 0.5 and 0.95 quantiles
Rscript -e 'quantile (as.numeric (readLines ("stdin")),probs = c(0.05, 0.1, 0.5, 0.90, 0.95))' < totalmissingsite
#anatina       1%        5%       50%       95%       99%
#           0.0012      0.025     0.092     0.440     0.488

awk '$1 <= 0.2' totalmissingsite | wc -l
awk '$1 <= 0.1' totalmissingsite | wc -l
# 52688864
# 38464153


# plot histogram
Rscript -e 'data <- as.numeric(readLines("stdin")); 
    pdf("aa_raw_miss5_Q30_DP2_bi_imiss5_site_miss.pdf"); 
    hist(data, main="Histogram of Missingness <0.5 per Site", xlab="Frequency missingness", col="steelblue", border="black", breaks=100); 
    dev.off()' < totalmissingsite

# ----------------------------

# ----------------------------
# ------- filt_2.sh ----------

VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5.vcf.gz"
# filter site based on missingess and mean site depth mindp 3 and max 0.99 quantiles and minor allele count of 2
vcftools --gzvcf ${VCF_FILE} --max-missing 0.8 --min-meanDP 3 --max-meanDP 13.568 --mac 2 --recode-INFO-all --recode --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss8_meanDP_mac2.vcf.gz}
# 34845398 SNPs

# filter for stranbias
bcftools filter -e 'INFO/PV4[0]< 0.01' ${VCF_FILE/.vcf.gz/_miss8_meanDP_mac2.vcf.gz} -Oz -o ${VCF_FILE/.vcf.gz/_miss8_meanDP_mac2_SB.vcf.gz}

# ----------------------------

# ----------------------------
# -------- stats.sh ----------

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


# ------- check strand bias ----------

VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5_miss9_meanDP_mac2.vcf.gz"

## STRAND BIAS ##

bcftools query -f '%PV4\n' ${VCF_FILE} > pv4.txt &
cut -d "," -f 1 pv4.txt > sb_pval &

awk '$1 >=0.01' sb_pval| wc -l 
# 29398903

awk '$1 >=0.05' sb_pval| wc -l 
# 26745690

# after bonferroni correction (PV4 < 0.01/34845398)
awk '$1 >=0.000000001435' sb_pval| wc -l
# 33665975

#plot histogram
Rscript -e 'data <- as.numeric(readLines("stdin")); 
    pdf("aa_raw_miss5_Q30_DP2_bi_imiss5_miss9_meanDP_mac2_sb_pval_hist.pdf"); 
    hist(data, main="Histogram of Strand Bias per Site", xlab="p-values", col="steelblue", border="black", breaks=100); 
    dev.off()' < sb_pval

# fdr correction
Rscript -e 'data <- as.numeric(readLines("stdin")); write.table(p.adjust(data, method = "BH"), "pvalues_fdr.txt", row.names = FALSE, col.names = "FDR_pval", quote = FALSE)'< sb_pval

awk '$1 >=0.01' pvalues_fdr.txt| wc -l 
# 31058179


# ------- Plot stats ----------

source $GDCstack
module load r

VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB.vcf.gz"
VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB_miss85_imiss3_mac5_renamed.vcf.gz"

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

# ----------------------------
# ------- filt_3.sh ----------

VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB.vcf.gz"

# filter for missingness
vcftools --gzvcf ${VCF_FILE} --max-missing 0.85  --recode-INFO-all --recode --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss85.vcf.gz} 

echo "Finished filtering and calculating individual stats $(date)"
vcftools --gzvcf ${VCF_FILE/.vcf.gz/_miss85.vcf.gz} --missing-indv --out ${VCF_FILE/.vcf.gz/_miss85} 

tail +2 ${VCF_FILE/.vcf.gz/_miss85}.imiss | awk '$5>0.3' | cut -f 1 > samples_imiss3

vcftools --gzvcf ${VCF_FILE/.vcf.gz/_miss85.vcf.gz} --remove samples_imiss3 --recode --recode-INFO-all --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss85_imiss3_mac5.vcf.gz} 


# vcf filters:
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

# Rename individual names in the vcf file
VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB_miss85_imiss3_mac5.vcf.gz"
# get old name and make file for converting
bcftools query -l $VCF_FILE > sample_names.txt
awk -F'[_-]' '{print $0, $3}' sample_names.txt > new_names.txt
# convert to uncompressed format
bcftools view -O v -o temp_uncompressed.vcf ${VCF_FILE}
# renames
bcftools reheader -s new_names.txt -o temp_renamed.vcf temp_uncompressed.vcf
# bgzip -c → Compresses the renamed VCF.
bgzip -c temp_renamed.vcf > ${VCF_FILE/.vcf.gz/_renamed.vcf.gz}
# index 
tabix -p vcf ${VCF_FILE/.vcf.gz/_renamed.vcf.gz}


# ----------------------------

# ----------------------------
# ------- subsampling ----------

# make 
# - bed files
# - smaller vcfs 
# - species subest

source $GDCstack
module load plink
module load vcftools/0.1.16-tc6l6nq
#mkdir beds

VCF_FILE="aa_raw_miss5_Q30_DP2_bi_imiss5_miss8_meanDP_mac2_SB_miss85_imiss3_mac5_renamed.vcf.gz"
BED_FILE="anatina"

#  make a bed file
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --out beds/${BED_FILE} &

# make a vcf file which has less info field and thus is smaller and easier to read into memory
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --recode vcf-iid --out beds/$BED_FILE &


awk -F',' '$1 == "Anodonta anatina" && $7 == "Ticino" {print $NF "\t" $NF}' metadata_all.csv  > anatina_TI_samples
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_TI_samples --make-bed --pca --out anatina_TI
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_TI_samples --mac 2 --make-bed --pca --out anatina_TI_mac2

awk -F',' '$1 == "Anodonta anatina" && $7 != "Ticino" {print $NF "\t" $NF}' metadata_all.csv  > anatina_N_samples
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_N_samples --make-bed --pca --out anatina_N
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_N_samples --mac 5 --make-bed --pca --out anatina_N_mac5 &
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_N_samples --recode vcf-iid --out anatina_N_mac5 &
plink --bfile $BED_FILE --allow-extra-chr --pca --out ${BED_FILE}

# ----------------------------

# ----------------------------
# -------- pcaone_ld.sh ----------

# filter for LD but account for population (species) structure

source $GDCstack
module load plink

# split vcf file in two for memory usage
BED_FILE="anatina"

ll $HOME/software/PCAone/ -b $BED_FILE -k 3 -n 1 --ld-r2 0.1 --ld-bp 100000 --out $BED_FILE

awk '{print $2 }' ${BED_FILE}.ld.prune.out > snp2rm

# Finally filter out SNPs in LD
#  make a bed file with ld filter ( here I also do pca )
plink --bfile $BED_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --exclude snp2rm --make-bed --pca --out ${BED_FILE}_ld

# make vcf file from filtered bed file
plink --bfile --out ${BED_FILE}_ld --allow-extra-chr --recode vcf-iid --out ${BED_FILE}_ld


BED_FILE="anatina_N_mac5"
# make vcf file from filtered bed file
plink --bfile $BED_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --exclude snp2rm --make-bed --pca --out ${BED_FILE}_ld

BED_FILE="anatina_TI_mac2"
# make vcf file from filtered bed file
plink --bfile $BED_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --exclude snp2rm --make-bed --pca --out ${BED_FILE}_ld



# ----------------------------





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
# quality.sh


# ----------------------------
# ------------AF-------------- DONE
# ----------------------------
out="afreq"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink2

BED_FILE="anatina"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out ${out}/$BED_FILE &

BED_FILE="anatina_N_mac5"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out ${out}/$BED_FILE &

BED_FILE="anatina_TI_mac2"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out ${out}/$BED_FILE &

# ----------------------------
# --------LD decay------------ DONE
# ----------------------------

out="ld"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink

BED_FILE="anatina"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.25 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --make-bed --out ${out}/${BED_FILE}_maf_thin
awk '{print $7, $5-$2}' ${out}/${BED_FILE}_maf_thin.ld > ${out}/${BED_FILE}_maf_thin_ld_with_distance.txt 

BED_FILE="anatina_N_mac5"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.25 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --make-bed --out ${out}/${BED_FILE}_maf_thin
awk '{print $7, $5-$2}' ${out}/${BED_FILE}_maf_thin.ld > ${out}/${BED_FILE}_maf_thin_ld_with_distance.txt 

BED_FILE="anatina_TI_mac2"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.25 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --make-bed --out ${out}/${BED_FILE}_maf_thin
awk '{print $7, $5-$2}' ${out}/${BED_FILE}_maf_thin.ld > ${out}/${BED_FILE}_maf_thin_ld_with_distance.txt 

# This command creates a files with
# r² value (R2)
# distance between POS_B and POS_A


# ==================================================================================================================================================
# 1 Population structure 
# ==================================================================================================================================================


# ----------------------------
# ------------PCA------------- DONE
# ----------------------------
# pca.sh

out="pca"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink

BED_FILE="anatina_ld"
#BED_FILE="anatina_N_mac5_ld"

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
# admix.sh

BED_FILE="anatina_ld"

cat ${BED_FILE}.bim | cut -f 2 > sites
shuf sites -n 2000000 | sort > sites2M
plink --bfile $BED_FILE --allow-extra-chr --extract sites2M --make-bed --out ${BED_FILE}_2M

BED_FILE="anatina_ld_2M"
BED_FILE="anatina_N_mac5_ld_2M"
# make folder and move files
source $GDCstack
module load plink


out="admix"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

# move to folder and copy in bed files
cd admix

cp -v ../${BED_FILE}.* ./





# ADMIXTURE does not accept chromosome names that are not human chromosomes. 
# We will thus just exchange the first column by 0

awk '{$1="0";print $0}' ${BED_FILE}.bim > ${BED_FILE}.bim.tmp
mv ${BED_FILE}.bim.tmp ${BED_FILE}.bim


# test run
#source $GDCstack
#module load admixture
#admixture --cv=5 -j8 ${BED_FILE}.bed 2 


BED_FILE="anatina_N_mac5_ld"
JOBID1=$(sbatch admix.sh $BED_FILE | awk '{print $4}')



BED_FILE="anatina_ld"
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


BED_FILE="anatina_ld"
BED_FILE="anatina"


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
plink2 --bfile ${BED_FILE} --double-id --allow-extra-chr --pheno pheno_5 --fst Pop --out ${out}/${BED_FILE}_pop5 &

# ----------------------------
# ---------RAxML----------
# ----------------------------

# https://github.com/amkozlov/raxml-ng/wiki/Tutorial 

BED_FILE="anatina_ld"
BED_FILE="anatina_ld_2M"

source $GDCstack
module load raxml-ng/1.2.2 
module load python
module load plink
module load bcftools

plink --bfile $BED_FILE --allow-extra-chr --recode vcf-iid --out ${BED_FILE}

#to leave only SNPs that were present at least once in our dataset as  homozygous for the reference allele, and homozygous for the alternative allele, as required by RAxML.
bcftools view -i 'COUNT(GT="RR")>0 & COUNT(GT="AA")>0' ${BED_FILE}.vcf -Ov -o ${BED_FILE}_filtered.vcf & # Uncompressed VCF
bcftools view -H ${BED_FILE}_filtered.vcf | wc -l 
# 1124070


#cat $BED_FILE.bim | cut -f 2 > sites
#shuf sites -n 100000 | sort > sites100000
#plink --bfile $BED_FILE --allow-extra-chr --extract sites100000 --recode vcf-iid --out ${BED_FILE}_100k


BED_FILE="anatina_ld_2M_filtered"
PHY_FILE="anatina_ld_2M_filtered.min4.phy"

# convert vcf format to phylip for our raxml analyses.
python $HOME/software/vcf2phylip-master/vcf2phylip.py -i ${BED_FILE}.vcf --output-folder ${out}


# check that the MSA can actually be read
raxml-ng --check --msa $PHY_FILE --model GTR+ASC_LEWIS --prefix T1

# Compress alignment patterns and store MSA in the binary format (RAxML Binary Alignment, RBA
# Estimate memory requirements and optimal number of CPUs/threads
raxml-ng --parse --msa $PHY_FILE --model GTR+ASC_LEWIS --prefix T2

raxml-ng --msa T2.raxml.rba --model GTR+ASC_LEWIS --prefix T3 --threads 22 --seed 2 # takes about 3mi











# Now let's infer a tree under GTR+GAMMA with default parameters. We will use 12 threads
# This command will perform 20 tree searches using 10 random and 10 parsimony-based starting trees, and pick the best-scoring topology:
raxml-ng --msa T2.raxml.rba --model GTR+G --prefix T3 --threads 12 --seed 2 # timed out after 24 h
#raxml-ng --all --msa T2.raxml.rba --model GTR+G --prefix T3 --threads 12 --seed 2 --tree pars{10} --bs-trees 100
# raxml-ng --msa oT2.raxml.rba --model GTR+G # for just one ML tree

# --all means do a full analysis: estimate best ML tree, then do bootstraps, then use bootstraps to plot support on best tree
# --msa means multiple-sequence alignment
# --model indicates evolutionary model, here we use General Time Reversible plus the Lewis ascertainment bias correction (takes into account the fact that we are only using variable sites)
# --bs-trees indicates the number of bootstraps you'd like to do


grep "Final LogLikelihood:" T3.raxml.log
# Final LogLikelihood: -1683690.144150
grep "ML tree search #" T3.raxml.log
# Here, we see why it is so important to use multiple starting trees: some searches converged to a local optimum with a much lower likelihood

raxml-ng --rfdist --tree T3.raxml.mlTrees --prefix RF3
# Average absolute RF distance in this tree set: 307.147368
# Average relative RF distance in this tree set: 0.386835
# Number of unique topologies in this tree set: 20

# ==================================================================================================================================================
# 2. Diversity
# ==================================================================================================================================================


# ----------------------------
# ------------dxy------------- 
# ----------------------------



# --------------------------------
# ---nucleotide diversity -------- DONE
# --------------------------------
# pi.sh


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

BED_FILE="anatina"



# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 && $9 == "passed" {print $NF "\t" $NF "\t" $2 "_" $5}' metadata_all.csv > pheno
# Extract unique groups (populations) from the third column of the pheno file
tail -n +2 pheno | cut -f 3 | sort | uniq > pops

grep "AA_" pops > pops_AA
# Individual heterozygosity to calculate FIS


while read -r group; do
    echo "working with group $group"
    awk -v grp="$group" '$3 == grp {print $1 "\t" $2 }' pheno > ${out}/keep$group
    #plink --bfile $BED_FILE --allow-extra-chr --keep keep$group --het --out ${BED_FILE}${group}
    vcftools --vcf ${BED_FILE}.vcf --keep ${out}/keep$group --het --out ${out}/${BED_FILE}${group}
done < pops_AA

# concatenate the files


echo -e "pop\tsample\tO(HOM)\tE(HOM)\tN_SITES\tF" > ${out}/${BED_FILE}_pop_het.txt

while read -r group; do
tail -n +2 ${out}/*${group}.het | awk -v grp="$group" '{print grp "\t" $0}' >> ${out}/${BED_FILE}_pop_het.txt
done < pops_AA


# ----------------------------
# ----------- ROH ------------ DONE
# ----------------------------

out="roh"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load bcftools
module load vcftools/0.1.16-tc6l6nq


VCF_FILE="anatina.vcf"

bgzip $VCF_FILE
tabix -p vcf ${VCF_FILE}.gz

# add allele frequency to file
bcftools +fill-tags ${VCF_FILE}.gz -Oz -o ${VCF_FILE/.vcf/}_AF.vcf.gz -- -t all &

# This limits ROH calculations to sites with 30 or more genotyped individuals.
# use the AF field for allele frequency
# only iutput RG
bcftools roh --GTs-only 30 --AF-tag AF ${VCF_FILE/.vcf/}_AF.vcf.gz --output-type r -o ${out}/${VCF_FILE}_roh.txt


# ----------------------------
# -------relatedness ---------
# ----------------------------
out="kinship"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

source $GDCstack
module load plink2


# List of BED files
BED_FILES=("anatina_N_mac5_ld" "anatina_TI_mac2_ld" "anatina_ld" "anatina" "anatina_N_mac5" "anatina_TI_mac2")

# Loop over BED files
for BED_FILE in "${BED_FILES[@]}"; do
    echo "Processing $BED_FILE..."
    
    # Run KING kinship estimator
    plink2 --bfile $BED_FILE --allow-extra-chr --make-king-table counts --make-king square --out ${out}/$BED_FILE
    
    # Check for high kinship values
    kin_file="${out}/${BED_FILE}.kin0"
    if [ -f "$kin_file" ]; then
        count=$(awk '$NF > 0.2' "$kin_file" | wc -l)
        echo "$kin_file: $count"
    else
        echo "Warning: $kin_file not found!"
    fi
done

# ==================================================================================================================================================
# 4. Effective population size
# ==================================================================================================================================================

# ----------------------------
# -------- currentNe --------- DONE
# ----------------------------

out="ne"
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi


VCF_FILE=anatina_N_mac5
# subset vcf file
bcftools view -H ${VCF_FILE}.vcf | shuf -n 1999999 | cut -f1,2 > ${out}/${VCF_FILE}_2Msnsp.txt

# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 && $9 == "passed" {print $NF "\t" $NF "\t" $2 "_" $5}' metadata_all.csv > pheno

# Extract unique groups (populations) from the third column of the pheno file
tail -n +2 pheno | cut -f 3 | sort | uniq -c | awk '$1 >= 10 {print $2}' > pops_10

grep "AA" pops_10 > pops_10_AA


# ne.sh

VCF_FILE=anatina_N_mac5
pops=pops_10_AA
out=ne


# how to reference each sample by array number
IDX=$SLURM_ARRAY_TASK_ID
group=`sed -n ${IDX}p < $pops`

echo "Processing file: $pops"

source $GDCstack
module load bcftools

# Create the 'keep' file for the group
awk -v grp="$group" '$3 == grp {print $1}' pheno > ${out}/keep$group
    
bcftools view -T ${out}/${VCF_FILE}_2Msnsp.txt -S ${out}/keep${group} -o ${out}/${group}.vcf -O v ${VCF_FILE}.vcf

echo "Staring Ne calculations of $pops"

$HOME/software/currentNe/currentNe -t 8 ${out}/${group}.vcf 19


# Once finished runs this on output

# To see Ne estimate an CI
grep -A9 "# Ne point estimate:" *_currentNe_OUTPUT.txt

# To only see Ne estimate
grep -A1 "# Ne point estimate:" *_currentNe_OUTPUT.txt

# Ne estimates without genetic map
for group in $(< pops_10_AA); do
    ne=$(awk 'NR==50' "${group}_currentNe_OUTPUT.txt")
    CI_10=$(awk 'NR==56' "${group}_currentNe_OUTPUT.txt")
    CI_90=$(awk 'NR==58' "${group}_currentNe_OUTPUT.txt")
    echo "$group $ne CI:${CI_10}-${CI_90}"
done


# AA_AGE 17.60 CI:14.14-21.91
# AA_ARH 11.60 CI:9.31-14.46
# AA_BRE 40.60 CI:32.36-50.94
# AA_BRU 23.95 CI:19.17-29.93
# AA_BUR 12.80 CI:9.94-16.49
# AA_LAC 34.92 CI:23.99-50.81
# AA_NOT 20.81 CI:16.47-26.30
# AA_SCN 31.92 CI:23.30-43.72
# AA_SCZ 5.85 CI:5.01-7.31
# AA_SEE 31.36 CI:21.85-44.99
# AA_TUR 57.31 CI:42.33-77.59
# AA_VAG 42.82 CI:35.28-51.98
# AA_WAG 33.77 CI:25.11-45.40
# AA_WAL 33.53 CI:24.85-45.24


# ==================================================================================================================================================
# TEST ZONE
# ==================================================================================================================================================







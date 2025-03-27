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
squeue --all # same as without flag on Euler?
squeue -u jconrads # same as above on Euler?

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
cut -f 1,2 $index | sort -k2 -n -r > sorted_scaffolds
# prepare list of bam files
ls bams/*bam > bam_list

# split sorted scaffolds according to size
head -n 500 sorted_scaffolds | cut -f 1 > largest_scaffolds_500
tail -n +501 sorted_scaffolds | cut -f 1 > smaller_scaffolds
# run bcfcalls.sh with largest 500 scaffolds
# run bcfcall_loop.sh with remaining scaffolds

# bcfcalls.sh will do the follwing
# 1. generate genotype likelihoods with bcftools mpileup
#		skipping indels, with a max read-depth of 30
# 2. call genotypes with bcftools call
#		using multiallelic-caller method, only outputting variant sites
sbatch bcfcalls.sh
sbatch bcfcalls_loop.sh


# --------  copy raw bcf files to project ----------

# ==================================================================================================================================================
# 6. Filtering
# ==================================================================================================================================================

# ----------------------------
# ------- concat.sh ----------

VCF_FILE="raw.vcf.gz"
# Concatenate bcf files
bcftools concat --file-list bcf_list -Oz -o $VCF_FILE

# missingness per individual
vcftools --gzvcf ${VCF_FILE} --missing-indv --out ${VCF_FILE}

# mean coverage across the sites.
vcftools --gzvcf ${VCF_FILE} --site-depth --out ${VCF_FILE}

# individual depth
vcftools --gzvcf ${VCF_FILE} --depth --out ${VCF_FILE}
# ---------------------------



# ----------------------------
# ------ interactive ---------
# get samples with high missingness

VCF_FILE="raw.vcf"

# plot individual missingness as a histogram
mawk '!/IN/' ${VCF_FILE/.vcf/}.imiss | cut -f5 > totalmissing

gnuplot << \EOF
set terminal dumb size 200, 30
set autoscale
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

# calculate average missingness
awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' ${VCF_FILE/.vcf/}.imiss
# 0.31629

# list samples which have more than 75% missingness
cat ${VCF_FILE/.vcf/}.imiss | awk '$5>0.75' | cut -f 1 > samples_imiss75
# Also remove 328261_05-AA05_S17 (duplicate)

# Number of samples to be removed
wc -l samples_imiss75

# ---------FILTERING----------
# Run separately for each (combined, anatina, others), example here for combined vcf file
# ------- filt_2.sh ----------

VCF_FILE="raw.vcf.gz"

# filter for site missingness, mapping quality, minimum depth, biallelic sites and remove samples with missingness above cutoff
vcftools --gzvcf $VCF_FILE --max-missing 0.5 --minQ 30 --minDP 2 --max-alleles 2 --min-alleles 2 --remove samples_imiss75 --recode --recode-INFO-all --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss5_Q30_DP2_bi_imiss75.vcf.gz}
#COMMAND OUTPUT
#Run Time = X seconds

# stats
VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75.vcf.gz"
# calculate missingness per individual
vcftools --gzvcf ${VCF_FILE} --missing-indv --out ${VCF_FILE/.vcf.gz/}

# mean coverage across the sites.
vcftools --gzvcf ${VCF_FILE} --site-depth --out ${VCF_FILE/.vcf.gz/}

# individual depth
vcftools --gzvcf ${VCF_FILE} --depth --out ${VCF_FILE/.vcf.gz/}
#Run Time = X seconds
# ---------------------------

# ----------------------------
# ------ interactive ---------
# get depth statistics
module load r

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75.vcf.gz"


# mean coverage across sites
cut -f3 ${VCF_FILE/.vcf.gz/}.ldepth.mean | tail +2 > meandepthpersite
# calculate quantiles
Rscript -e 'quantile (as.numeric (readLines ("stdin")), probs = c(0.01, 0.05, 0.5, 0.95, 0.99))' < meandepthpersite
#anatina       1%        5%       50%       95%       99%
#        0.916667  1.333330  3.583330  6.291670 10.750000

#others       1%       5%      50%      95%      99%
#        1.77987  2.03774  4.55346 12.03770 29.92450

#combined      1%       5%      50%      95%      99%
#         1.03357  1.46043  3.81775  7.20624 13.29980


# calculate mean
Rscript -e 'mean (as.numeric (readLines ("stdin")))' < meandepthpersite
#anatina 4.775306
#others

# calculate sd
Rscript -e 'sd (as.numeric (readLines ("stdin")))' < meandepthpersite
#anatina
#2.013125
#2.013125 *3 = 6.039375
#others

#Let's plot the mean coverage across the sites.
gnuplot << \EOF
set terminal dumb size 200, 30
set autoscale
set xrange [0:25]
unset label
set title "Histogram of mean depth per site"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
set xtics 5
plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

#Let's check the missingness across the sites.
tail +2 ${VCF_FILE/.vcf.gz/}.lmiss | cut -f6 > totalmissingsite
# calculate 0.5 and 0.95 quantiles
Rscript -e 'quantile (as.numeric (readLines ("stdin")),probs = c(0.01, 0.05, 0.5, 0.95, 0.99))' < totalmissingsite
#anatina       1%        5%       50%       95%       99%
#           0.0012      0.025     0.092     0.440     0.488
#others       1%        5%       50%       95%       99%
#      0.00000000 0.00628931 0.10691800 0.42767300 0.48427700


# try 20-25 depth or even three standard deviation from the mean as max?

#Let's plot the missingness across the sites.
gnuplot << \EOF
set terminal dumb size 200, 30
set autoscale
set xrange [0:0.5]
unset label
set title "Histogram of missingness per site"
set ylabel "Number of Occurrences"
set xlabel "fraction missing per site"
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissingsite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

gnuplot << \EOF
set terminal dumb size 200, 30
set autoscale
#set xrange [0:0.5]
unset label
set title "Histogram of Strand bias Phred scale"
set ylabel "Number of Occurrences"
set xlabel "Phred scale"
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'sb_pval' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

# ----------------------------

# ----------------------------
# ------- filt_3.sh ----------

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75.vcf.gz"

# filter site missingess and mean site depth between 0.05 and 0.95 quantiles
vcftools --gzvcf ${VCF_FILE} --max-missing 0.9 --min-meanDP 2.32 --max-meanDP 7.22 --recode-INFO-all --recode --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_miss9_meanDP.vcf.gz}
# ----------------------------

# ----------------------------
# ------- filt_4.sh ----------

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP.vcf.gz"

# remove sites with strand bias p -value below 0.01 after bonferroni correction (PV4 < 0.01/74194034)
bcftools filter -e 'INFO/PV4[0]< 0.000000000134' $VCF_FILE -Oz -o ${VCF_FILE/.vcf.gz/_SB.vcf.gz}
# ----------------------------

# ----------------------------
# ------- filt_5.sh ----------

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB.vcf.gz"

# filter for minor allele count of 2
vcftools --gzvcf $VCF_FILE --mac 2 --recode --recode-INFO-all --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_mac2.vcf.gz}
# ----------------------------

# ----------------------------
# -------- stats.sh ----------

module load vcftools

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac2.vcf.gz"


echo "Finished filtering and calculating individual stats $(date)"
vcftools --gzvcf $VCF_FILE --depth --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --missing-indv --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --het --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --relatedness --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --relatedness2 --out ${VCF_FILE/.vcf.gz/}

echo "Calculating allele stats $(date)"
vcftools --gzvcf $VCF_FILE --freq2 --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --counts2 --out ${VCF_FILE/.vcf.gz/}

echo "Calculating site stats $(date)"
vcftools --gzvcf $VCF_FILE --site-depth --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --site-mean-depth --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --missing-site --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --singletons --out ${VCF_FILE/.vcf.gz/}
vcftools --gzvcf $VCF_FILE --hardy --out ${VCF_FILE/.vcf.gz/}


tail -n+2 ${VCF_FILE/.vcf.gz/.idepth} | awk -F '\t' '$3<3' | wc -l
# 416 individuals have an average depth below 10
# 149 individuals have an average depth below 5
# 33 individuals have an average depth below 3
# 3 individuals have an average depth below 1

tail -n+2 ${VCF_FILE/.vcf.gz/.imiss} | awk -F '\t' '$5>0.3' | wc -l
# 143 individuals with missingnes >5% 
# 72 individuals with missingnes >10%
# 47 individuals with missingnes >15%
# 35 individuals with missingnes >20%
# 27 individuals with missingnes >25%
# 17 individuals with missingnes >30%
# 11 individuals with missingnes >40%
# 8 individuals with missingnes >50%

tail -n+2 ${VCF_FILE/.vcf.gz/.lmiss}  | awk -F '\t' '$6>0.05' | wc -l
# 2 629 659 SNPs have a missingness above 5% (>21 individuals)
# 0 SNPs have a missingness above 10% (>42 individuals)

tail -n+2 ${VCF_FILE/.vcf.gz/.singletons} | wc -l
# 3469 SNPs only occur in one individual
tail -n+2 ${VCF_FILE/.vcf.gz/.singletons} | grep "D" | wc -l
# 3469 Doubletons

tail -n+2 ${VCF_FILE/.vcf.gz/.frq} | awk -F '\t' '$6<0.01' | wc -l
# 12 614 534 138 SNPs have a MAF below 0.05 (<42 occurences)
# 4 871 493 SNPs have a MAF below 0.01 (<9 occurences)


tail -n+2 ${VCF_FILE/.vcf.gz/.frq.count} | awk -F '\t' '$6<5' | wc -l
# 2 774 978 SNPs have a mac < 5



# ----------------------------
# ------- filt_6.sh ----------

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac2.vcf.gz"

tail -n+2 ${VCF_FILE/.vcf.gz/.imiss} | awk -F '\t' '$5>0.3'  > ind2rm

# filter for minor allele frequency and individual with >30 missing
vcftools --gzvcf $VCF_FILE --maf 0.01 --remove ind2rm --recode --recode-INFO-all --stdout | gzip -c > ${VCF_FILE/.vcf.gz/_mac2.vcf.gz}

# vcf filters:
# variant >50% missing data
# minQ 30
# individual with to much missing data  >0.75
# min DP 2
# bi-allelic
# variant <90% missing data
# mean DP between combined:1.46-7.21 
# strand bias p-value combined< 0.000000000305
# imiss >30%
# mac 5


# ------------------------------
# ------- subsampling ----------

# make 
# - bed files
# - smaller vcfs 
# - species subest

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load plink
module load vcftools/0.1.16-tc6l6nq

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac5_imiss30.vcf.gz"
BED_FILE="all_species"

#  make a bed file
plink2 --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --out ${BED_FILE}

# make a vcf file which has less info field and thus is smaller and easier to read into memory
plink2 --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --recode vcf-iid --out $BED_FILE


###### ANATAINA #########

# to create a tab delimited file for plink of individuals to keep

awk -F',' '$1 == "Anodonta anatina" {print $NF "\t" $NF}' groups_aa.csv  > anatina_samples
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_samples --make-bed --out anatina

awk -F',' '$1 == "Anodonta anatina" && $6 == "Ticino" {print $NF "\t" $NF}' groups_aa.csv  > anatina_TI_samples
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_TI_samples --make-bed --out anatina_TI
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_TI_samples --mac 2 --make-bed --recode vcf-iid --out anatina_TI_mac2 &

grep -vf anatina_TI_samples anatina_samples > anatina_N_samples
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_N_samples --make-bed --out anatina_N
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_N_samples --mac 5 --make-bed --recode vcf-iid --out anatina_N_mac5 &
###### CYGNEA #########

awk -F',' '$1 == "Anodonta cygnea" {print $NF "\t" $NF}' groups_aa.csv  > cygnea_samples
plink --bfile $BED_FILE --allow-extra-chr --keep cygnea_samples --make-bed --out cygnea

grep -vE "AC83|AC90" cygnea_samples > cygnea_samples_exc_hy
plink --bfile $BED_FILE --allow-extra-chr --keep cygnea_samples_exc_hy --mac 5 --make-bed --recode vcf-iid --out cygnea_mac5 &


###### EXULCERATA #########

awk -F',' '$1 == "Anodonta exulcerata" {print $NF "\t" $NF}' groups_aa.csv  > exulcerata_samples
plink --bfile $BED_FILE --allow-extra-chr --keep exulcerata_samples --make-bed --out exulcerata
plink --bfile $BED_FILE --allow-extra-chr --keep exulcerata_samples --mac 2 --make-bed --recode vcf-iid  --out exulcerata_mac2 &

###### EXULCERATA & CYGNEA #########

awk -F',' '$1 == "Anodonta exulcerata" {print $NF "\t" $NF}' groups_aa.csv  > exulcerata_samples
plink --bfile $BED_FILE --allow-extra-chr --remove anatina_samples --make-bed --out exulcerata_cygnea
plink --bfile $BED_FILE --allow-extra-chr --remove anatina_samples --recode vcf-iid  --out exulcerata_cygnea

# more subest options
awk -F',' '$1 == "Anodonta cygnea" && $5 == "UFS" {print $NF}' groups_aa.csv  > cygnea_UFS
awk -F',' '$1 == "Anodonta anatina" && $5 == "VAG" {print $NF}' groups_aa.csv  > anatina_VAG
# ----------------------------
# -------- pcaone_ld.sh ----------

# filter for LD but account for population (species) structure

module load stack
module load plink

# split vcf file in two for memory usage
BED_FILE="all_species"

./PCAone -b $BED_FILE -k 3 -n 1 --ld-r2 0.1 --ld-bp 100000 --out pcaone_all_species

awk '{print $2 }' pcaone_all_species.ld.prune.out > snp2rm

# Finally filter out SNPs in LD
#  make a bed file with ld filter ( here I also do pca )
plink --bfile $BED_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --exclude snp2rm --make-bed --out ${BED_FILE}_ld

# make vcf file from filtered bed file
plink --bfile --out ${BED_FILE}_ld --allow-extra-chr --recode vcf-iid --out ${BED_FILE}_ld


# ------------------------------
# ------- sub-sampling after LD ----------

# make 
# - bed files
# - smaller vcfs 
# - species subest

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load plink
module load vcftools/0.1.16-tc6l6nq
module load bcftools


BED_FILE="all_species_ld"


###### ANATAINA #########

# to create a tab delimited file for plink of individuals to keep
awk -F',' '$1 == "Anodonta anatina" {print $NF "\t" $NF}' groups_aa.csv  > anatina_samples
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_samples --mac 1 --pca --make-bed --out anatina_ld
plink --bfile $BED_FILE --allow-extra-chr --keep anatina_samples --mac 1 --recode vcf-iid  --out anatina_ld


###### CYGNEA #########

awk -F',' '$1 == "Anodonta cygnea" {print $NF "\t" $NF}' groups_aa.csv  > cygnea_samples
plink --bfile $BED_FILE --allow-extra-chr --keep cygnea_samples --mac 1 --pca --make-bed --out cygnea_ld
plink --bfile $BED_FILE --allow-extra-chr --keep cygnea_samples --mac 1 --recode vcf-iid  --out cygnea_ld

###### EXULCERATA #########

awk -F',' '$1 == "Anodonta exulcerata" {print $NF "\t" $NF}' groups_aa.csv  > exulcerata_samples
plink --bfile $BED_FILE --allow-extra-chr --keep exulcerata_samples --mac 1 --pca --make-bed --out exulcerata_ld
plink --bfile $BED_FILE --allow-extra-chr --keep exulcerata_samples --mac 1 --recode vcf-iid  --out exulcerata_ld

###### EXULCERATA & CYGNEA #########

awk -F',' '$1 == "Anodonta exulcerata" {print $NF "\t" $NF}' groups_aa.csv  > exulcerata_samples
plink --bfile $BED_FILE --allow-extra-chr --remove anatina_samples --mac 1 --pca --make-bed --out exulcerata_cygnea_ld
plink --bfile $BED_FILE --allow-extra-chr --remove anatina_samples --mac 1 --recode vcf-iid  --out exulcerata_cygnea_ld




# ------------------------------
# --------largest scaffolds-------------
# ------------------------------

source /cluster/project/gdc/shared/stack/GDCstack.sh

module load plink
module load vcftools/0.1.16-tc6l6nq
module load bcftools

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac5_imiss30.vcf.gz"

REF=./Aanat_demux_hifi_l3.bp.p_ctg.fa

# index
samtools faidx ${REF}

# Sort the .fai file by scaffold length and get the top 10:
sort -k2,2nr ${REF}.fai > sorted_scaffolds

awk '{sum += $2} END {print sum}' sorted_scaffolds
# 2896208859

head -10 sorted_scaffolds > top10_scaffolds # 0.6% of tot
head -20 sorted_scaffolds > top20_scaffolds # 1.1% of tot

head -20 sorted_scaffolds | cut -f1 > top20_scaffolds 
head -10 sorted_scaffolds | cut -f1 > top10_scaffolds 

# for the use of bcftools do:
#grep -f top20_scaffolds ${REF}.fai | awk '{print $1, ,1, ,$2}'  OFS="\t" > top20_scaffolds_inorder

vcftools --gzvcf $VCF_FILE $(awk '{print "--chr", $1}' top20_scaffolds | tr '\n' ' ') --recode --out all_species_top20
# After filtering, kept 253479 out of a possible 26704704 Sites

# to create a tab delimited file for plink of individuals to keep
awk -F',' '$1 == "Anodonta anatina" {print $NF}' groups_aa.csv  > anatina_samples

VCF_FILE="all_species_top20.recode.vcf"
vcftools --vcf $VCF_FILE --keep anatina_samples --maf 0.01 --recode --out anatina_top20
# After filtering, kept 245 out of 400 Individuals
# Outputting VCF file...
# After filtering, kept 91097 out of a possible 253479 Sites


vcftools --vcf $VCF_FILE --remove anatina_samples --maf 0.01 --recode --out cygnea_exulcerata_top20
# After filtering, kept 155 out of 400 Individuals
# Outputting VCF file...
# After filtering, kept 109152 out of a possible 253479 Sites


# ------------------------------
# --------phasing-------------
# ------------------------------

# cygnea & exulcerata
bgzip cygnea_exulcerata_top20.recode.vcf
VCF_FILE=cygnea_exulcerata_top20.recode.vcf.gz
tabix -p vcf $VCF_FILE

bcftools +fill-tags $VCF_FILE -Oz -o cygnea_exulcerata_top20_info.vcf.gz -- -t all

VCF_FILE="cygnea_exulcerata_top20_info.vcf.gz"
tabix -p vcf $VCF_FILE

$HOME/software/phase_common_static --input $VCF_FILE --region ptg001722l --output AC_ptg001722l.bcf


# anatina
bgzip anatina_top20.recode.vcf
VCF_FILE=anatina_top20.recode.vcf.gz
tabix -p vcf $VCF_FILE

bcftools +fill-tags $VCF_FILE -Oz -o anatina_top20_info.vcf.gz -- -t all

VCF_FILE="anatina_top20_info.vcf.gz"
tabix -p vcf $VCF_FILE

$HOME/software/phase_common_static --input $VCF_FILE --region ptg001722l --output AA_ptg001722l.bcf


# for for loop of phasing multiple regions see scripts
# phase.sh

# ------------------------------
# ------LD on haplotypes---------
# ------------------------------


BCF_FILE="AC_top20_phased.bcf"
BCF_FILE="AA_top20_phased.bcf"

BCF_FILE="AC_top20_phased.bcf"
vcftools --bcf $BCF_FILE --hap-r2 --ld-window-bp 500000 --out ${BCF_FILE/.bcf/} &
awk '{print $5, $3-$2}' ${BCF_FILE/.bcf/}.hap.ld > ${BCF_FILE/.bcf/}_ld_with_distance.txt &






awk '{print $5, $3-$2}' ${BCF_FILE/.bcf/}.hap.ld | head






# rehader to better names
grep "^#CHROM" combined.vcf | cut -f10- | tr '\t' '\n' | sed 's/^[^_]*_[^_]*-\([^_]*\)_.*/\1/' > newname\n
bcftools reheader --samples newname --output anadonta.vcf combined.vcf



#find specific files based of patterns
for zip in *.zip; do unzip -l "$zip" | grep '\.ldepth'; done

# to extract one specific file
unzip -j combined.zip combined/raw_miss5_Q30_DP2_bi_imiss75.ldepth.mean -d $SCRATCH/stats

# copy out specific files based of patterns
for zip in *.zip; do
    unzip -l "$zip" | grep '\.imiss' | awk '{print $NF}' | while read file; do
        unzip -j "$zip" "$file" -d $SCRATCH/stats
    done
done

unzip -l cygnea.zip | grep 'filt' | awk '{print $NF}' | while read file; do
        unzip -j cygnea.zip "$file" -d $SCRATCH/anodonta/cygnea/
    done

unzip -j others.zip others/raw_other_miss5_Q30_DP2_bi_imiss6_miss8_meanDP_SB_mac2.ldepth.mean -d $SCRATCH/anodonta/cygnea/



########### TEST ZONE #################





source /cluster/apps/local/env2lmod.sh
module load bcftools

out=/cluster/scratch/jconrads/SNPs
Ref=/cluster/scratch/jconrads/reference/Aanat_demux_hifi_l3.bp.p_ctg.fa
bams=/cluster/scratch/jconrads/var_call_job/bam_list

# if there is no out directory make one
if [ ! -e ${out} ]  ; then mkdir ${out} ; fi

#down-sample region with > max. cov. per individual 
maxcov=50
bcftools mpileup -f ${Ref} --skip-indels -R positions.txt -d ${maxcov} -a 'FORMAT/DP,FORMAT/AD' -b bam_list | bcftools call -mv -a -Ob  -o ${out}/raw.${name}.bcf

bcftools mpileup -f ${Ref} --skip-indels -R positions.txt -d ${maxcov} -a 'FORMAT/DP,FORMAT/AD' -b bam_list | bcftools call -mv -a 'INFO/PV4' -Ob  -o ${out}/raw.${name}.bcf


head -n 200 raw_allele_depth_aa.txt | grep --color=always -E "0/1[:]4[,]0" | head -n1


awk -F',' '$1 == "Anodonta anatina" {print $NF "_sort_Q20_nodup.bam"}' groups_aa.csv  > anatina_bams



cat raw_allele_depth_aa.txt | grep --color=always -E "0/1[:]0[,][0-9]+|0/1[:][0-9][,]0+" | head -n10000 | cut -f1 | sed 's/:/\t/g' >aa_bad_pos.txt






### OLD STUFF



#----------------------------
# -------- LD filter --------

# Create a linkage pruned data set
module load plink
module load vcftools
module load bcftools

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac2.vcf.gz"
BED_FILE="combined"

#  make a bed file
#plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --out ${BED_FILE}_23

# perform linkage pruning - i.e. identify prune sites
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out linkage
# Pruning complete.  18674708 of 21626081 variants removed. (2.951373 million SNPs left)

#  make a bed file with ld filter ( here I also do pca )
plink --vcf $VCF_FILE --double-id --allow-extra-chr --set-missing-var-ids @:# --extract linkage.prune.in --make-bed --pca --out $BED_FILE

# make vcf file from filtered bed file
plink --bfile $BED_FILE --allow-extra-chr --recode vcf-iid --out $BED_FILE


# ----------------------------
# -------- 10k filter --------

# make a down sampled data set for quick diversity estimates

# make random samples list, did this for 10k, 100k and 1M
cat anatina_noprune.bim | cut -f 2 > sites
shuf sites -n 10000 | sort > sites10000
plink --bfile $BED_FILE --allow-extra-chr --extract sites10000 --recode vcf-iid --out ${BED_FILE/2.9/10k}


# if you want to keep all the data from vcf file you can also
# create snp list for filtering with vcf
#sed s'/:/\t/'g linkage.prune.in > snps2keep

# creating ld filtered vcf file
#vcftools --gzvcf $VCF_FILE --positions snps2keep --recode --recode-INFO-all --out raw_miss5_Q30_DP3_bi_imiss5_SB_miss9_meanDP_norep_mac_ld.vcf



# to calculate per cluster

vcftools --gzvcf $VCF_FILE --site-pi --out ${VCF_FILE/.vcf.gz/} # per site pi
vcftools --gzvcf $VCF_FILE --hardy --out ${VCF_FILE/.vcf.gz/} # observed and expected H per site and p-values
vcftools --gzvcf $VCF_FILE --SNPdensity 1000 --out ${VCF_FILE/.vcf.gz/} # density of snps in bins on 1kb
vcftools --gzvcf $VCF_FILE --counts2 --out ${VCF_FILE/.vcf.gz/} # get counts, can look for monomorphic sites and number of polymorphic sites

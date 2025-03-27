# Population analysis pipeline for duckmussel Anodonta anatina
# Last updated by Ellika Faust September 2024

# ==================================================================================================================================================
# Data analysis
# ==================================================================================================================================================

# 0. Quality
# 1. Population structure (PCA, admixture, Fst)
# 2. Genetic diversity (Ho, He, Polymorphic sites, D, Pi, Theta)
# 3. Inbreeding (Fis, IBS, ROH, genetic load)
# 4. Effective population size (Ne, PSMC)



# T chain jobs:
JOBID1=$(sbatch job1.sh | awk '{print $4}')
JOBID2=$(sbatch --dependency=afterany:$JOBID1 admix.sh | awk '{print $4}')
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
# ------------AF-------------
# ----------------------------
source /cluster/project/gdc/shared/stack/GDCstack.sh
module load plink2


BED_FILE="anatina"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out $BED_FILE

BED_FILE="anatina_N"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out $BED_FILE

BED_FILE="anatina_TI"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out $BED_FILE

BED_FILE="cygnea"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out $BED_FILE

BED_FILE="exulcerata"
plink2 --bfile $BED_FILE --allow-extra-chr --freq --out $BED_FILE




# ----------------------------
# ------------LD decay-------------
# ----------------------------
module load stack/2024-06  gcc/12.2.0
module load plink

BED_FILE="exulcerata_pol"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.1 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --out $BED_FILE
awk '{print $7, $5-$2}' ${BED_FILE}.ld > ${BED_FILE}_ld_with_distance.txt &

BED_FILE="cygnea_pol"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.1 --r2 --ld-window 100 --ld-window-kb 100  --ld-window-r2 0 --out $BED_FILE 
awk '{print $7, $5-$2}' ${BED_FILE}.ld > ${BED_FILE}_ld_with_distance.txt &

BED_FILE="anatina_pol"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.1 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --out $BED_FILE
awk '{print $7, $5-$2}' ${BED_FILE}.ld > ${BED_FILE}_ld_with_distance.txt &

BED_FILE="anatina_TI_pol"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.1 --r2 --ld-window 100 --ld-window-kb 100  --ld-window-r2 0 --out $BED_FILE
awk '{print $7, $5-$2}' ${BED_FILE}.ld > ${BED_FILE}_ld_with_distance.txt &

BED_FILE="anatina_N_pol"
plink --bfile $BED_FILE --allow-extra-chr --maf 0.01 --thin 0.1 --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --out $BED_FILE
awk '{print $7, $5-$2}' ${BED_FILE}.ld > ${BED_FILE}_ld_with_distance.txt &


# This command creates a files with
# r² value (R2)
# distance between POS_B and POS_A


# ----------------------------
# ------------DP per GT-------------
# ----------------------------
source /cluster/project/gdc/shared/stack/GDCstack.sh
module load bcftools

VCF_FILE="all_species_top20.recode.vcf"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%DP]\n' $VCF_FILE > depth_per_sample.tsv

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT,%DP]\n' $VCF_FILE > genotype_depth.tsv

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%AD]\n' $VCF_FILE | head

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%DP]\n' $VCF_FILE| head -n 2

VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac5_imiss30_anatina_mac5.recode.vcf"

VCF_FILE="cygnea_exulcerata_top20_info.vcf.gz"
VCF_FILE="anatina_top20_info.vcf.gz"
bcftools query -f '%CHROM:%POS\t[%GT:%DP;]\n' $VCF_FILE > raw_genotype_depth.txt


awk '
BEGIN {
    print "SNP\tn_00\tn_01\tn_11\tDP_00\tDP_01\tDP_11"
}
{
    split($2, samples, ";");
    sum_00 = sum_01 = sum_11 = 0;
    count_00 = count_01 = count_11 = 0;
    for (i in samples) {
        split(samples[i], info, ":");
        genotype = info[1];
        depth = info[2];
        if (genotype == "0/0") {
            sum_00 += depth;
            count_00++;
        } else if (genotype == "0/1") {
            sum_01 += depth;
            count_01++;
        } else if (genotype == "1/1") {
            sum_11 += depth;
            count_11++;
        }
    }
    avg_00 = (count_00 > 0) ? sum_00 / count_00 : 0;
    avg_01 = (count_01 > 0) ? sum_01 / count_01 : 0;
    avg_11 = (count_11 > 0) ? sum_11 / count_11 : 0;
    print $1 "\t" avg_00 "\t" avg_01 "\t" avg_11 "\t" count_00 "\t" count_01 "\t" count_11;
}' raw_genotype_depth.txt > snp_depth_avg_table.txt


VCF_FILE="cygnea_exulcerata_top20_info.vcf.gz"
VCF_FILE="anatina_top20_info.vcf.gz"
bcftools query -f '%CHROM:%POS\t[%GT:%AD;]\n' $VCF_FILE > raw_genotype_allele_depth.txt





awk '
BEGIN {
    print "SNP\tn_01\tDP_0\tDP_1"
}
{
    split($2, samples, ";");
    sum_01 = 0;
    count_01 = 0;
    sum_01_allele1 = 0;
    sum_01_allele0 = 0;
    
    for (i in samples) {
        split(samples[i], info, ":");
        genotype = info[1];
        allele_depth = info[2];
        split(allele_depth, depths, ",");
        depth_0 = depths[1];  # Depth for allele 0
        depth_1 = depths[2];  # Depth for allele 1

        # Handling genotypes
        if (genotype == "0/1") {
            sum_01_allele0 += depth_0;
            sum_01_allele1 += depth_1;
            count_01++;
        }
    }

    # Print results for each SNP
    print $1 "\t" count_01 "\t" sum_01_allele0 "\t" sum_01_allele1;
}' raw_genotype_allele_depth.txt > cygnea_ad.txt




awk '
BEGIN {
    print "SNP\tn_01\tDP_0\tDP_1"
}
{
    split($2, samples, ";");
    sum_01_allele0 = 0;
    sum_01_allele1 = 0;
    count_01 = 0;

    for (i in samples) {
        split(samples[i], info, ":");
        genotype = info[1];

        # Skip missing genotypes and homozygotes
        if (genotype == "./." || genotype == "0/0" || genotype == "1/1") {
            continue;
        }
        allele_depth = info[2];
        split(allele_depth, depths, ",");
        depth_0 = depths[1];  # Depth for allele 0
        depth_1 = depths[2];  # Depth for allele 1

        # Process heterozygotes only
        if (genotype == "0/1") {
            sum_01_allele0 += depth_0;
            sum_01_allele1 += depth_1;
            count_01++;
        }
    }

    # Print results for each SNP
    print $1 "\t" count_01 "\t" sum_01_allele0 "\t" sum_01_allele1;
}' raw_allele_depth_ac.txt > cygnea_ad.txt



bcftools filter -e 'GT!="0/1" || AD[0] < 1 || AD[1] < 1' -Oz -o filtered.vcf.gz input.vcf.gz


# ==================================================================================================================================================
# 1 Population structure 
# ==================================================================================================================================================


# ----------------------------
# ------------PCA-------------
# ----------------------------


source /cluster/project/gdc/shared/stack/GDCstack.sh
module load plink

BED_FILE="all_species_ld"





# ----------------------------
# ---------ADMIXTURE----------
# ----------------------------


BED_FILE="all_species_ld"
BED_FILE="anatina_ld"
BED_FILE="cygnea_ld"
BED_FILE="exulcerata_cygnea_ld"


# make folder and move files
mkdir admix
cp -v ${BED_FILE}.* ./admix/
cd admix

# ADMIXTURE does not accept chromosome names that are not human chromosomes. 
# We will thus just exchange the first column by 0

awk '{$1="0";print $0}' ${BED_FILE}.bim > ${BED_FILE}.bim.tmp
mv ${BED_FILE}.bim.tmp ${BED_FILE}.bim


# test run
#source /cluster/project/gdc/shared/stack/GDCstack.sh
#module load admixture
#admixture --cv=5 -j8 ${BED_FILE}.bed 2 


# move bed files into admixture folder and run
sbatch admix.sh

# note that admixture needs a lot of memory for high Ks, especially when calculating cv for many SNPs
# if it keeps running out of memory try thinning the data, either by LD or by position.

# get cross validation values

awk ' /CV/ {print $3,$4}' *${BED_FILE}*out | sed -e 's/(K=//;s/)://' | sort -n -k 1 > ${BED_FILE}.cv





# ----------------------------
# ---------    FST  ----------
# ----------------------------



source /cluster/project/gdc/shared/stack/GDCstack.sh
module load plink2
BED_FILE="all_species_ld"
BED_FILE="cygnea_mac5"



# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 {print $NF "\t" $NF "\t" $2 "_" $5}' groups_aa.csv > pheno

plink2 --bfile ${BED_FILE} --double-id --allow-extra-chr --pheno pheno --fst Pop --out ${BED_FILE}_pop



# ----------------------------
# ---------RAxML----------
# ----------------------------

# https://github.com/amkozlov/raxml-ng/wiki/Tutorial 


source /cluster/project/gdc/shared/stack/GDCstack.sh
module load raxml-ng/1.2.2 
module load python

BED_FILE="all_species"
cat $BED_FILE.bim | cut -f 2 > sites
shuf sites -n 100000 | sort > sites100000
plink --bfile $BED_FILE --allow-extra-chr --extract sites100000 --recode vcf-iid --out ${BED_FILE}_100k


# convert vcf format to phylip for our raxml analyses.
python $HOME/software/vcf2phylip-master/vcf2phylip.py -i all_species_100k.vcf


# check that the MSA can actually be read
raxml-ng --check --msa out.phy --model GTR+G --prefix T1

# Compress alignment patterns and store MSA in the binary format (RAxML Binary Alignment, RBA
# Estimate memory requirements and optimal number of CPUs/threads
raxml-ng --parse --msa out.phy --model GTR+G --prefix T2

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

# this is a python script to calculate dxy between to files of allele frequencies 

cat << EOF > calculate_average_dxy.py
import pandas as pd
import argparse

def calculate_average_dxy(pop1_file, pop2_file):
    # Load population data
    pop1 = pd.read_csv(pop1_file, sep="\t")
    pop2 = pd.read_csv(pop2_file, sep="\t")

    # Merge the two populations on the ID column
    merged = pd.merge(pop1, pop2, on="ID", suffixes=("_pop1", "_pop2"))

    # Calculate dxy for each site
    merged['dxy'] = (
        merged['ALT_FREQS_pop1'] * (1 - merged['ALT_FREQS_pop2']) +
        (1 - merged['ALT_FREQS_pop1']) * merged['ALT_FREQS_pop2']
    )

    # Calculate the average dxy across all sites
    average_dxy = merged['dxy'].mean()

    # Output the result
    print(f"Average dxy: {average_dxy}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate the average dxy between two populations")
    parser.add_argument("pop1_file", help="File for population 1")
    parser.add_argument("pop2_file", help="File for population 2")

    # Parse command line arguments
    args = parser.parse_args()

    # Run the average dxy calculation
    calculate_average_dxy(args.pop1_file, args.pop2_file)
EOF


source /cluster/project/gdc/shared/stack/GDCstack.sh
module load python

# initialize output files
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

# Add a header for the new column and process the file
awk -v frac="$FRACTION" '
    NR == 1 {print $0 "\tCorrected_dxy"; next}  # Add header for the new column
    NR > 1 {print $0 "\t" $NF * frac}          # Process data rows
' all_dxy_results.txt > updated_all_dxy_results.txt

# ----------------------------
# ---------nucleotide diversity -------------
# ----------------------------

# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 {print $NF "\t" $NF "\t" $2 "_" $5}' groups_aa.csv > pheno
# Extract unique groups (populations) from the third column of the pheno file
tail -n +2 pheno | cut -f 3 | sort | uniq > pops

# Individual heterozygosity to calculate FIS
VCF_FILE="raw_miss5_Q30_DP2_bi_imiss75_miss9_meanDP_SB_mac5_imiss30.vcf.gz"

while read -r group; do
    # Create the 'keep' file for the group
    awk -v grp="$group" '$3 == grp {print $1}' pheno > keep$group
    
    # Check the number of lines (individuals) in the 'keep' file
    num_individuals=$(wc -l < keep$group)

    # Only run vcftools if there are more than 5 individuals
    if [ "$num_individuals" -gt 4 ]; then
        echo "Running vcftools for group $group with $num_individuals individuals"
        vcftools --gzvcf $VCF_FILE --window-pi 10000 --keep keep$group --out pi_10kb_$group
    else
        echo "Skipping group $group, not enough individuals ($num_individuals)"
    fi
done < pops


# concatenate the files
head -n1 pi_50kb_AE_AGN.windowed.pi | awk -v grp="GROUP" '{print grp "\t" $0}' > group_pi.txt
while read -r group; do
tail -n +2 pi_50kb_${group}.windowed.pi | awk -v grp="$group" '{print grp "\t" $0}' >> group_pi.txt
done < pops


# ----------------------------
# -----------Het------------
# ----------------------------

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load plink
module load vcftools/0.1.16-tc6l6nq


grep "AE" pops > pops_AE
grep "AA" pops | grep -v "AA_LOC" > pops_AA
grep "AC" pops > pops_AC


BED_FILE="exulcerata_mac2"

while read -r group; do
    echo "working with group $group"
    awk -v grp="$group" '$3 == grp {print $1 "\t" $2 }' pheno > keep$group
    plink --bfile $BED_FILE --allow-extra-chr --keep keep$group --het --out ${BED_FILE}${group}
    vcftools --vcf ${BED_FILE}.vcf --keep keep$group --het --out vcftools_${BED_FILE}${group}
done < pops_AE

BED_FILE="cygnea_mac5"

while read -r group; do
    echo "working with group $group"
    awk -v grp="$group" '$3 == grp {print $1 "\t" $2 }' pheno > keep$group
    plink --bfile $BED_FILE --allow-extra-chr --keep keep$group --het --out ${BED_FILE}${group}
    vcftools --vcf ${BED_FILE}.vcf --keep keep$group --het --out vcftools_${BED_FILE}${group}
done < pops_AC

BED_FILE="anatina_N_mac5"

while read -r group; do
    echo "working with group $group"
    awk -v grp="$group" '$3 == grp {print $1 "\t" $2 }' pheno > keep$group
    plink --bfile $BED_FILE --allow-extra-chr --keep keep$group --het --out ${BED_FILE}${group}
    vcftools --vcf ${BED_FILE}.vcf --keep keep$group --het --out vcftools_${BED_FILE}${group}
done < pops_AA

BED_FILE="anatina_TI_mac2"
group=AA_LOC
plink --bfile $BED_FILE --allow-extra-chr --het --out ${BED_FILE}${group}
vcftools --vcf ${BED_FILE}.vcf --het --out vcftools_${BED_FILE}${group}


# concatenate the files
head -n1 cygnea_mac5AC_WAD.het | awk -v grp="pop" '{print grp "\t" $0}' > sep_species_pop_het.txt
while read -r group; do
tail -n +2 *${group}.het | awk -v grp="$group" '{print grp "\t" $0}' >> sep_species_pop_het.txt
done < pops

head -n1 cygnea_mac5.het | awk -v grp="pop" '{print grp "\t" $0}' > sep_species_het.txt
# Loop over each species
for group in cygnea_mac5 exulcerata_mac2 anatina_N_mac5 anatina_TI_mac2; do
    # Append the data with the group name as the first column
    tail -n +2 "${group}.het" | awk -v grp="$group" '{print grp "\t" $0}' >> sep_species_het.txt
done

# ==================================================================================================================================================
# 3. Inbreeding
# ==================================================================================================================================================


# ----------------------------
# -----------FIS------------
# ----------------------------



# make pheno file with pop info
awk -F, 'BEGIN {OFS="\t"; print "FID" "\t" "IID" "\t" "Pop"} NR > 1 {print $NF "\t" $NF "\t" $2 "_" $5}' groups_aa.csv > pheno
# Extract unique groups (populations) from the third column of the pheno file
tail -n +2 pheno | cut -f 3 | sort | uniq > pops

# Individual heterozygosity to calculate FIS
BED_FILE="all_species"


while read -r group; do
    echo "working with group $group"
    awk -v grp="$group" '$3 == grp {print $1 "\t" $2 }' pheno > keep$group
    plink --bfile $BED_FILE --allow-extra-chr --keep keep$group --het --out ${BED_FILE}${group}
    vcftools --vcf ${BED_FILE}.vcf --keep keep$group --het --out vcftools_${BED_FILE}${group}
done < pops_AC


plink --bfile $BED_FILE --allow-extra-chr --het --out ${BED_FILE}
vcftools --vcf ${BED_FILE}.vcf --het --out vcftools_${BED_FILE}



# ----------------------------
# ----------- ROH ------------
# ----------------------------

# bcftools

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load bcftools
module load vcftools/0.1.16-tc6l6nq


VCF_FILE="all_species_top20.recode.vcf"
bgzip $VCF_FILE
tabix -p vcf ${VCF_FILE}.gz
bcftools +fill-tags ${VCF_FILE}.gz -Oz -o all_species_top20_info.vcf.gz -- -t all
VCF_FILE="all_species"
bcftools roh --GTs-only 30 ${VCF_FILE}_top20_info.vcf.gz -o ${VCF_FILE}_top20_roh.txt 




awk -F',' '$1 == "Anodonta anatina" {print $NF "\t" $NF}' groups_aa.csv  > anatina_samples
awk -F',' '$1 == "Anodonta anatina" && $6 == "Ticino" {print $NF "\t" $NF}' groups_aa.csv  > anatina_TI_samples
grep -vf anatina_TI_samples anatina_samples > anatina_N_samples
awk -F',' '$1 == "Anodonta cygnea" {print $NF "\t" $NF}' groups_aa.csv  > cygnea_samples
awk -F',' '$1 == "Anodonta exulcerata" {print $NF "\t" $NF}' groups_aa.csv  > exulcerata_samples


for VCF_FILE in all_species anatina_N anatina_TI cygnea exulcerata; do
    VCF_FILE_ALL="all_species_top20.recode.vcf.gz"
    #vcftools --gzvcf $VCF_FILE_ALL --keep ${VCF_FILE}_samples --recode --out ${VCF_FILE}_top20
    #bgzip ${VCF_FILE}_top20.recode.vcf
    #tabix -p vcf ${VCF_FILE}_top20.recode.vcf.gz
    #bcftools +fill-tags ${VCF_FILE}_top20.recode.vcf.gz -Oz -o ${VCF_FILE}_top20_info.vcf.gz -- -t all
    bcftools roh --GTs-only 30 ${VCF_FILE}_top20_info.vcf.gz -O r -o ${VCF_FILE}_top20_roh.txt 
done


echo -e "RG\tSample\tChromosome\tStart\tEnd\tLength\tN\tQuality" > merged_top20_roh.txt
for VCF_FILE in anatina_N anatina_TI cygnea exulcerata; do
    grep "^RG" ${VCF_FILE}_top20_roh.txt >> merged_top20_roh.txt
done

echo -e "RG\tSample\tChromosome\tStart\tEnd\tLength\tN\tQuality" > all_species_top20_roh_clean.txt
grep "^RG"  all_species_top20_roh.txt >> all_species_top20_roh_clean.txt




# ----------------------------
# -----------relatedness ------------
# ----------------------------


source /cluster/project/gdc/shared/stack/GDCstack.sh
module load plink2
BED_FILE="all_species_ld"
BED_FILE="cygnea_ld"
BED_FILE="anatina_ld"
BED_FILE="exulcerata_ld"
BED_FILE="exulcerata_cygnea_ld"


# KING kinship estimator
plink2 --bfile $BED_FILE --allow-extra-chr --make-king-table counts --make-king square --out $BED_FILE

# to check  for high values
awk '$NF > 0.25' *.kin0 | wc -l



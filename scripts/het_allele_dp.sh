## Get allele depth for heterozygots 

VCF_FILE="anatina_top20_info.vcf.gz"
bcftools query -f '%CHROM:%POS\t[%GT:%AD;]\n' $VCF_FILE > raw_genotype_allele_depth.txt

# Iterates Over Each SNP and Sample:
## Skips missing genotypes (./.) and homozygous genotypes (0/0, 1/1).
## Extracts allele depths (AD) for heterozygous genotypes (0/1).
## Adds the sequencing depth of allele 0 to sum_01_allele0.
## Adds the sequencing depth of allele 1 to sum_01_allele1.
## Increments the heterozygous sample count (count_01).

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


# Plot in R

ad_all <- read_table("cygnea_ad.txt") %>%
  filter(n_01!=0) %>% # filter out sites with 0 het calls
  mutate(f_1=DP_0/n_01, # get average depth
         f_0= DP_1/n_01,
         diff=abs(f_1-f_0)) # get absolut difference between alles

ad_all %>% 
  ggplot(aes(x=diff)) + 
  geom_density(alpha = .4) +
  labs(x = "difference between average allele depth across all heterozygotes", y ="")

ggsave(filename = "het_allele_dp.png")
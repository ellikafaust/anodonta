#!/bin/bash
#SBATCH --job-name=phase
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00
#SBATCH --output=%x.o%A
#SBATCH --error=%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch


###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$


# Path to the input VCF file
VCF_FILE="anatina_top20_info.vcf.gz"

# Path to the regions file
REGIONS_FILE="top20_scaffolds"

# Path to the software
PHASE_SOFTWARE="$HOME/software/phase_common_static"

# Loop through each region listed in the file
while read -r region; do
    echo "Processing region: $region"
    
    # Run the command for each region
    $PHASE_SOFTWARE --input $VCF_FILE --region "$region" --output "AA_${region}.bcf"
    
    echo "Finished processing region: $region"
done < "$REGIONS_FILE"

ls -1v AA*bcf > AA_bcf_files
bcftools concat --naive -f AA_bcf_files -o AA_top20_phased.bcf 




#########################################
#----------------------------------------
#########################################

# Path to the input VCF file
VCF_FILE="cygnea_exulcerata_top20_info.vcf.gz"

# Loop through each region listed in the file
while read -r region; do
    echo "Processing region: $region"
    
    # Run the command for each region
    $PHASE_SOFTWARE --input $VCF_FILE --region "$region" --output "AC_${region}.bcf"
    
    echo "Finished processing region: $region"
done < "$REGIONS_FILE"

ls -1v AC*bcf > bcf_files
bcftools concat -f bcf_files -o AC_top20_phased.bcf 


#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$


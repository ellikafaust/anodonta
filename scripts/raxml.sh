#!/bin/bash
#SBATCH --job-name=raxml
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
#SBATCH --output=%x.o%A
#SBATCH --error=%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch


###############################################################################################################################################$
echo "Starting ${SLURM_ARRAY_TASK_ID} at $(date)"
###############################################################################################################################################$

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load raxml-ng
module load python

# convert vcf format to phylip for our raxml analyses.
#python $HOME/software/vcf2phylip-master/vcf2phylip.py -i all_species_100k.vcf
 
raxml-ng --msa T2.raxml.rba --model GTR+ASC_LEWIS --prefix T3 --threads auto{12} --seed 2 # takes about 3mi


#################################
echo "Job: ${SLURM_ARRAY_TASK_ID} successfully finished $(date)"
# happy end
exit 0
###############################################################################################################################################$

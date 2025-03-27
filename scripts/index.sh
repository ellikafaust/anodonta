#!/bin/bash
#SBATCH --job-name=indexv3 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G
#SBATCH --time=1:00:00
#SBATCH --output=%x.o%A
#SBATCH --error=%x.e%A
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ellika.faust@eawag.ch

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load bwa-mem2
bwa-mem2 index ./ref/v3.asm.bp.p_ctg.fa.gz
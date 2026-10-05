#!/bin/bash
#SBATCH -A hpc2ncourses2026-019
#SBATCH -J map_guides_hg19
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=map_guides_hg19_%j.out
#SBATCH --error=map_guides_hg19_%j.err

module load GCC/13.3.0
module load Bowtie/1.3.1

bowtie -f -v 0 -a -p 1 hg19 library_b_targeting_guides.fa guide_hg19_exact.tsv

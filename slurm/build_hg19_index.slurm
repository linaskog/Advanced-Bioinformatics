#!/bin/bash
#SBATCH -A hpc2ncourses2026-019
#SBATCH -J build_hg19_index
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=build_hg19_index_%j.out
#SBATCH --error=build_hg19_index_%j.err

module load GCC/13.3.0
module load Bowtie/1.3.1

bowtie-build hg19.fa hg19


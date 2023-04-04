#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=run_rseqc
#SBATCH --time=08:00:00
#SBATCH --mem=8G

mamba activate rseqc

geneBody_coverage.py -i bampaths.txt -r GRCm39_chr1-4.bed -o rseqc_qc/group;

exit

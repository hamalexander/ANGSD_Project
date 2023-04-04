#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=index_genome
#SBATCH --time=05:00:00
#SBATCH --mem=60G

mamba activate angsd

STAR --runMode genomeGenerate \
--runThreadN 1 \
--genomeDir GRCm39_STARindex \
--genomeFastaFiles GRCm39_genome.fa \
--sjdbGTFfile GRCm39.gtf \
--sjdbOverhang 99

exit

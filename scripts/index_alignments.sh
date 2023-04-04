#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=index_alignments
#SBATCH --time=02:00:00
#SBATCH --mem=8G

mamba activate angsd

for dir in ERR*/; do
        filename=$(find ${dir} -maxdepth 1 -name "*.bam" -print -quit | xargs basename);
        path=${dir}${filename};
	samtools index ${path};
done

exit

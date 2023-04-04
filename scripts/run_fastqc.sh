#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=run_fastqc
#SBATCH --time=02:00:00
#SBATCH --mem=8G

mamba activate angsd

for file in *.fastq.gz; do
	fastqc ${file} -o fastqc_reports/;
done

exit
	

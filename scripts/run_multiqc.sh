#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=run_multiqc
#SBATCH --time=02:00:00
#SBATCH --mem=8G

mamba activate multiqc

for dir in ERR*/; do
	filename=$(cut -d "_" -f 3 <<< ${dir})
	multiqc --outdir ./multiqc_reports --filename "${filename}_multiqc" ./${dir};
done

exit

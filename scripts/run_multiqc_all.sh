#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=run_multiqc_all
#SBATCH --time=02:00:00
#SBATCH --mem=8G

mamba activate multiqc

dirs=(ERR*/)

multiqc --outdir ./multiqc_reports --filename "all_multiqc" "${dirs[@]}"

exit

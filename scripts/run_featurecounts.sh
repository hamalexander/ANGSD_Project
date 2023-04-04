#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=run_featurecounts
#SBATCH --time=02:00:00
#SBATCH --mem=8G

mamba activate angsd

path_array=()

for dir in ERR*/; do
        filename=$(find ${dir} -maxdepth 1 -name "*.bam" -print -quit | xargs basename);
        path=${dir}${filename};
        if [[ "${dir}" != "ERR3094492_Beta_LFD1/" ]]; then
		path_array+=($path);
	fi;
done

featureCounts -a GRCm39.gtf \
-o project_featureCounts.txt \
-t exon \
-g gene_id \
-p \
--countReadPairs \
-O \
-T 4 \
${path_array[@]};

exit


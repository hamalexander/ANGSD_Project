#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=align_reads
#SBATCH --time=02:00:00
#SBATCH --mem=36G

mamba activate angsd

pair1="None"
pair2="None" 
for file in *.fastq.gz; do
	if [[ ${pair1} == "None" ]]; then 
		pair1=${file}
	else
		pair2=${file}
		dir_name=$(cut -d "_" -f 1-3 <<< ${pair1}) 
		mkdir ${dir_name}

		STAR --runMode alignReads \
		--runThreadN 8 \
		--genomeDir GRCm39_STARindex \
		--readFilesIn ${pair1} ${pair2} \
		--readFilesCommand zcat \
		--outFileNamePrefix ./${dir_name}/${dir_name} \
		--outSAMtype BAM SortedByCoordinate \

		pair1="None"
		pair2="None"
	fi;
done
exit
	 

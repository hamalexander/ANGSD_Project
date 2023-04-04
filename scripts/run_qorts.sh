#!/bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=run_qorts
#SBATCH --time=02:00:00
#SBATCH --mem=8G

spack find | grep -i qorts
spack load qorts@1.2.42
QORTS_LOC=$(spack location -i qorts)

for dir in ERR*/; do
        filename=$(find ${dir} -maxdepth 1 -name "*.bam" -print -quit | xargs basename);
        path=${dir}${filename};
        
	java -Xmx8G -jar ${QORTS_LOC}/bin/QoRTs.jar QC \
	--generatePdfReport \
	${path} \
	GRCm39_v2.gtf.gz \
	qorts_qc
done

exit

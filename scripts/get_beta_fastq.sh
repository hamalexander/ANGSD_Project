#!/bin/bash
#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=get_beta_fastq
#SBATCH --time=02:00:00
#SBATCH --mem=4G 

while read line; do
	run=$(echo "$line" | cut -f2)
	link=$(echo "$line" | cut -f3)
 	link1="ftp://"$(echo "$link" | cut -d';' -f1)
        link2="ftp://"$(echo "$link" | cut -d';' -f2)
	condition=$(echo ${line} | egrep "Beta_...." -o)
	
	filename=${run}"_"${condition}
	echo ${link1} | xargs wget -O ${filename}_1.fastq.gz 
	echo ${link2} | xargs wget -O ${filename}_2.fastq.gz

 
done < beta_links.tsv

exit

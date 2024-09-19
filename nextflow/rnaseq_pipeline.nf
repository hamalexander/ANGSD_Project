#!/usr/bin/env nextflow

params.outdir = "./results"
params.genomeDir = "GRCm39_STARindex"
params.gtf = "GRCm39.gtf"

process fastqc {
    container 'rnaseq-pipeline'
    
    input:
    file fastq from Channel.fromFilePairs("*.fastq.gz")

    output:
    file "fastqc_reports/*" into fastqc_reports

    script:
    """
    fastqc ${fastq} -o fastqc_reports/
    """
}

process alignReads {
    container 'rnaseq-pipeline'

    input:
    file fastqPairs from Channel.fromFilePairs("*.fastq.gz")
    path genomeDir from params.genomeDir

    output:
    file "*.bam" into bamFiles

    script:
    """
    pair1=\$(echo ${fastqPairs[0]})
    pair2=\$(echo ${fastqPairs[1]})

    dir_name=\$(basename \$pair1 .fastq.gz)
    mkdir -p \$dir_name

    STAR --runMode alignReads \
         --runThreadN 8 \
         --genomeDir ${genomeDir} \
         --readFilesIn \$pair1 \$pair2 \
         --readFilesCommand zcat \
         --outFileNamePrefix ./${dir_name}/${dir_name} \
         --outSAMtype BAM SortedByCoordinate
    """
}

process runMultiQC {
    container 'rnaseq-pipeline'

    input:
    file bamDirs from bamFiles

    output:
    file "multiqc_reports/*" into multiqcReports

    script:
    """
    multiqc --outdir ./multiqc_reports --filename "all_multiqc" "${bamDirs}"
    """
}

process runFeatureCounts {
    container 'rnaseq-pipeline'

    input:
    file bamFiles from bamFiles
    path gtf from params.gtf

    output:
    file "project_featureCounts.txt" into featureCountsOutput

    script:
    """
    featureCounts -a ${gtf} \
                  -o project_featureCounts.txt \
                  -t exon \
                  -g gene_id \
                  -p \
                  --countReadPairs \
                  -O \
                  -T 4 \
                  ${bamFiles}
    """
}

workflow {
    fastqc()
    alignReads()
    runMultiQC()
    runFeatureCounts()
}

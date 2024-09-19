#!/usr/bin/env nextflow

params.genomeDir = "GRCm39_STARindex"
params.genomeFasta = "GRCm39_genome.fa"
params.gtf = "GRCm39.gtf"

process indexGenome {
    container 'rnaseq-pipeline'

    input:
    path genomeFasta from params.genomeFasta
    path gtf from params.gtf

    output:
    path "${params.genomeDir}" into genomeIndex

    script:
    """
    STAR --runMode genomeGenerate \
         --runThreadN 1 \
         --genomeDir ${params.genomeDir} \
         --genomeFastaFiles ${genomeFasta} \
         --sjdbGTFfile ${gtf} \
         --sjdbOverhang 99
    """
}

workflow {
    indexGenome()
}
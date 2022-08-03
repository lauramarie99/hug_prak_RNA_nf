#! /usr/bin/env nextflow
// Define parameters

params.reads = "~/praktikum/data/gut_{1,2}.fq"
params.transcriptome_file = "~/praktikum/data/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "$projectDir/results"

log.info """\
    NEXTFLOW RNA-SEQ PIPELINE
    =========================
    transcriptome   :   ${params.transcriptome_file}
    reads           :   ${params.reads}
    outdir          :   ${params.outdir}
    """
    .stripIndent()

println "ProjectDir: ${projectDir}"

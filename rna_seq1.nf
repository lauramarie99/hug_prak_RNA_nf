#! /usr/bin/env nextflow
# comment
params.reads = "$projectDir/data/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"

println "reads: $params.reads"

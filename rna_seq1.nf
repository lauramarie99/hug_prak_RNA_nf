#! /usr/bin/env nextflow

// Define parameters
params.reads = '/home/ec2-user/praktikum/data/*_{1,3}.fq'
params.transcriptome_file = "/home/ec2-user/praktikum/data/transcriptome.fa"
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

// Create binary index of transcriptome file
process index {
    cpus 2
    conda "/home/ec2-user/anaconda3/envs/salmon_env"
    
    input:
    path transcriptome from params.transcriptome_file

    output:
    path "salmon_index" into index_ch

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .set {read_pairs_ch}

read_pairs_ch.view()


#! /usr/bin/env nextflow

// Define parameters
params.reads = '/home/ec2-user/praktikum/data/*_{1,2}.fq'
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

// Read Pairs Channel
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    //.set {read_pairs_ch}
    .into {read_pairs_ch; read_pairs_ch2}


// Quantification
process quantification {
    conda "/home/ec2-user/anaconda3/envs/salmon_env"
    tag "$pair_id"
    publishDir "${params.outdir}/quantification", mode: "copy"

    input:
    path index from index_ch
    tuple pair_id, path(reads) from read_pairs_ch

    output:
    path pair_id into quant_ch

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

// Quality control using FastQC
process fastqc {
    conda "/home/ec2-user/anaconda3/envs/fastqc"
    tag "fastqc on ${pair_id}"

    input:
    tuple pair_id, path(reads) from read_pairs_ch2

    output:
    path "fastqc_${pair_id}_logs" into fastqc_ch

    script:
    """
    mkdir fastqc_${pair_id}_logs
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${reads}
    """
}

// Multi QC 
process multiqc {

    conda "/home/ec2-user/anaconda3/envs/multiqc"
    publishDir params.outdir, mode: 'copy'

    input:
    path '*' from quant_ch.mix(fastqc_ch).collect()

    output:
    path 'mulitqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow.onComplete {
    log.info(workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops.. something went wrong")
}

// pipeline input parameters
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.reads = "$projectDir/data/samplesheet.csv"
params.outdir = "results"

/*
 * define the `INDEX` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {

    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
    publishDir params.outdir, mode: 'copy'

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index -t $transcriptome -i salmon_index
    """
}

process QUANTIFICATION {

    tag "salmon on ${sample_id}"
    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path salmon_index
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --libType=U -i $salmon_index -1 ${reads_1} -2 ${reads_2} -o $sample_id
    """
}

process FASTQC {

    tag "fastqc on ${sample_id}"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir -p "fastqc_${sample_id}_logs"
    fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2 -t $task.cpus
    """
}

process MULTIQC {

    container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
    publishDir params.outdir, mode: 'copy'

    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
    fastqc_ch = FASTQC(read_pairs_ch)
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)

    fastqc_ch
        .mix(quant_ch)
        .collect()
        .set { all_qc_ch }

    MULTIQC(all_qc_ch)
}

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

process MULTIQC {

    container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
    publishDir params.outdir, mode: 'copy'

    input:
    path "*"

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}

// Define the workflow
workflow {

  // Run the index step with the transcriptome parameter
  INDEX(params.transcriptome_file)

  // Define the fastqc input channel
  reads_in = Channel.fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }

  // Run the fastqc step with the reads_in channel
  FASTQC(reads_in)

  // Run the quantification step with the index and reads_in channels
  transcriptome_index_in = INDEX.out[0]
  QUANTIFICATION(transcriptome_index_in, reads_in)

  // Define the multiqc input channel
  multiqc_in = FASTQC.out[0]
    .mix(QUANTIFICATION.out[0])
    .collect()

  /*
   * Generate the analysis report with the 
   * outputs from FASTQC and QUANTIFICATION
   */ 
  MULTIQC(multiqc_in)

}

/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.reads = "$projectDir/data/samplesheet.csv"
params.outdir = "results"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `INDEX` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
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

workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .set{ read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(read_pairs_ch)
}

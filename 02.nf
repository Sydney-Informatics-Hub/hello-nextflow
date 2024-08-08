/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.outdir = "results"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)

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

workflow {
    Channel
        .fromFilePairs(params.reads)
        .set { read_pairs_ch }
    read_pairs_ch.view()
    index_ch = INDEX(params.transcriptome_file)
}

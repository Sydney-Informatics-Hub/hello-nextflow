/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
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

workflow {
    index_ch = INDEX(params.transcriptome_file)
}

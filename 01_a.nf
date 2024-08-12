/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"

/*
 * define the `INDEX` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {

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

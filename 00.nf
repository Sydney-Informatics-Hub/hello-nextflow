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

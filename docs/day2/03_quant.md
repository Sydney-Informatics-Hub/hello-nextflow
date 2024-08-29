# 2.2 Perform expression quantification

## 2.2.1 Implementing the process  

Add the following `QUANTIFICATION` process definition to the script:  

```groovy title="main.nf"
process QUANTIFICATION {

    input:
    path salmon_index
    tuple val(sample_id), path(reads_1), path(reads_2)
    
    output:
    path "$sample_id"
    
    script:
    """
    salmon quant --libType=U -i $salmon_index -1 ${reads_1} -2 ${reads_2} -o $sample_id
    """
```

!!! question "Exercise"

    Add the `container` and `publishDir` directives to the `QUANTIFICATION` process.  
    ??? note "Solution"

        ```groovy title="main.nf"
        process QUANTIFICATION {

            container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
            publishDir params.outdir, mode: 'copy'
            
            input:
        ```

The process requires 4 inputs: 
- `salmon_index` - the output of the `INDEX` process  
- `sample_id` - the sample name/prefix of the paired `.fq` reads  
- `reads_1/2` - the relative file paths to the paired `.fq` reads  

You have already successfully ran `salmon` and obtained the required `output`
and now will prepare the remaining inputs for the process.  

Lastly, we will assign the channel output to a variable using the
[`set`](https://www.nextflow.io/docs/latest/operator.html#set) operator. `set`
can be used to define a new channel variable in place of an `=` assignment.  

Remove `.view()` and assign the channel output with `set`:  

```groovy title="main.nf"
workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)]
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
```

Update the workflow with the following:  

```groovy title="main.nf"
    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(read_pairs_ch)
}
```

Run the workflow:  
```bash
nextflow run main.nf
```

```console title="Output"
Launching `main.nf` [shrivelled_cuvier] DSL2 - revision: 4781bf6c41

executor >  local (2)
[35/05252c] INDEX              | 1 of 1 ✔
[dc/b8763a] QUANTIFICATION (1) | 1 of 1 ✔

```

!!! question "Exercise"

    What is the name of the directory that was output by `QUANTIFICATION`?

    ??? note "Solution"

    `gut` or `results/gut`

> Mention output, cache, resume  

!!! abstract "Summary"

    In this step you have learned:

        1. How to          
        1. How to 
        1. How to 
        1. How to 
        1. How to 

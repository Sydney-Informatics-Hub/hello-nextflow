# 2.2 Perform expression quantification

## 2.2.1 Read file pairs with a samplesheet  

https://nextflow-io.github.io/patterns/process-per-csv-record/

Inspect samplesheet:  

```bash
cat data/samplesheet.csv
```

```console, title="Output"
sample, fastq_1, fastq_2
gut, data/ggal/gut_1.fq, data/ggal/gut_2.fq
```

The samplesheet consists of three columns:  

- `sample`: indicates the sample name/prefix  
- `fastq_1`, `fastq_2` contains the relative paths to the reads  

!!! question "Exercise"  

    Add a parameter called `reads` and assign the path to the samplesheet.  

    ??? note "Solution"  

        ```groovy title="main.nf"
        /*
         * pipeline input parameters
         */
        params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
        params.reads = "$projectDir/data/samplesheet.csv"
        params.outdir = "results"
        ```

!!! question "Exercise"

    Add `params.reads` to `log.info`.  

    ??? Solution

        ```groovy title="main.nf"
        log.info """\
            R N A S E Q - N F   P I P E L I N E
            ===================================
            transcriptome: ${params.transcriptome_file}
            reads        : ${params.reads}
            outdir       : ${params.outdir}
            """
            .stripIndent()
        ```

Use some groovy to read and view the samplesheet. Define the channel in the
workflow, and add the
[`view`](https://www.nextflow.io/docs/latest/operator.html#view) operator to see
the contents of the channel:  

```groovy title="main.nf"
workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .view()

    index_ch = INDEX(params.transcriptome_file)
```

Add graphic.  

[`Channel.fromPath`](https://www.nextflow.io/docs/latest/channel.html#frompath)
creates a channel from the file path (i.e. reads in the samplesheet).

The [splitCsv](https://www.nextflow.io/docs/latest/operator.html#splitcsv)
operator splits the samplesheet line-by-line. The `header: true` argument
indicates that the first row is a header.  

Run the workflow.

```bash
nextflow run main.nf  
```

Your output will display a tuple consisting of three named elements, one for
each of the columns in the samplesheet.  

```console title="Output"
Launching `main.nf` [ecstatic_goldwasser] DSL2 - revision: 2de2dd7b2a

R N A S E Q - N F   P I P E L I N E
===================================
transcriptome: /home/ubuntu/hello-nextflow/data/ggal/transcriptome.fa
reads        : /home/ubuntu/hello-nextflow/data/samplesheet.csv
outdir       : results

executor >  local (1)
[20/0727ea] INDEX | 1 of 1 ✔
[sample:gut,  fastq_1: data/ggal/gut_1.fq,  fastq_2: data/ggal/gut_2.fq]

```

The first element of the tuple is the sample name (`gut`), and `fastq_1` and
`fastq_2` are the paths to the read files.  

!!! note

    We will explore how this works for multiple samples (rows in the
    samplesheet) in the later steps.  

Next, we will assign the channel output to a variable using the
[`set`](https://www.nextflow.io/docs/latest/operator.html#set) operator. `set`
can be used to define a new channel variable in place of an `=` assignment.  

Remove `.view()` and assign the channel output with `set`:  

```groovy title="main.nf"
workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
```

## 2.2.2 Implementing the process  

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

    Add the `publishDir` directive to the `QUANTIFICATION` process.  

    ??? note "Solution"

        ```groovy title="main.nf"
        process QUANTIFICATION {
            publishDir params.outdir, mode: 'copy'
            
            input:
        ```

Update the workflow with the following:  

```groovy title="main.nf"
workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(params.transcriptome_file)
}
```

Run the workflow:  
```bash
nextflow run main.nf
```

!!! question "Poll"

    1. What was the error recieved?  
    2. See the docs on [multiple input channels](https://www.nextflow.io/docs/latest/process.html#multiple-input-channels). Which scope should be amended to fix the issue?
        a. `process { }`
        b. `workflow { }`



Add `tag` for a more readable execution log. Will also help with profiling
later, when additional samples are added.  

> Question: Why would a `tag` directive not be added to `INDEX`?  

Run the workflow.
```
Launching `main.nf` [reverent_lavoisier] DSL2 - revision: 10860f201c

R N A S E Q - N F   P I P E L I N E
===================================
transcriptome: /home/fredjaya/GitHub/hello-nextflow/data/ggal/transcriptome.fa
reads        : /home/fredjaya/GitHub/hello-nextflow/data/ggal/gut_{1,2}.fq
outdir       : results

executor >  local (2)
[11/010b59] INDEX                          | 1 of 1 ✔
[21/5e9ec8] QUANTIFICATION (salmon on gut) | 1 of 1 ✔
```

> Question: Check `results/`, what is the name of the output directory?  

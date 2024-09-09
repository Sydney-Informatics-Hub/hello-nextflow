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

## 2.2.1 Read file pairs with a samplesheet  

> Prose about operators  

> Prose about groovy

> Prose about sample sheets  

https://nextflow-io.github.io/patterns/process-per-csv-record/

Inspect samplesheet:  

```bash
cat data/samplesheet.csv
```

```console, title="Output"
sample,fastq_1,fastq_2
gut,data/ggal/gut_1.fq,data/ggal/gut_2.fq
```

The samplesheet consists of three columns:  

- `sample`: indicates the sample name/prefix  
- `fastq_1`, `fastq_2` contains the relative paths to the reads  

!!! question "Exercise"  

    Add a parameter called `reads` and assign the path to the samplesheet using
    `$projectDir`.  

    ??? note "Solution"  

        ```groovy title="main.nf"
        /*
         * pipeline input parameters
         */
        params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
        params.reads = "$projectDir/data/samplesheet.csv"
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

> Add graphic  

[`Channel.fromPath`](https://www.nextflow.io/docs/latest/channel.html#frompath)
creates a channel from the file path (i.e. reads in the samplesheet).

The [splitCsv](https://www.nextflow.io/docs/latest/operator.html#splitcsv)
operator splits the samplesheet line-by-line. The `header: true` argument
indicates that the first row is a header.  

Before you execute the workflow again, define `outdir` in the script to avoid
calling it in the command line each time.  

!!! question "Exercise"

    Assign `"results"` to `params.outdir` in the workflow script.

    ??? note "Solution"  

        ```groovy title="main.nf"
        /*
         * pipeline input parameters
         */
        params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
        params.reads = "$projectDir/data/samplesheet.csv"
        params.outdir = "results"
        ```

Run the workflow.

```bash
nextflow run main.nf  
```

Your output will display a tuple consisting of three named elements, one for
each of the columns in the samplesheet.  

```console title="Output"
Launching `main.nf` [ecstatic_goldwasser] DSL2 - revision: 2de2dd7b2a

executor >  local (1)
[20/0727ea] INDEX | 1 of 1 ✔
[sample:gut, fastq_1:data/ggal/gut_1.fq, fastq_2:data/ggal/gut_2.fq]

```

The first element of the tuple is the sample name (`gut`), and `fastq_1` and
`fastq_2` are the paths to the read files.  

!!! note

    We will explore how this works for multiple samples (rows in the
    samplesheet) in the later steps.  

One more command needs to be added to ensure that the inputs are in the correct
format for the process. Add the following:  

```groovy title="main.nf"
workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)]
        .view()

    index_ch = INDEX(params.transcriptome_file)
```

Run the workflow.

```bash
nextflow run main.nf
```

Your output should look similar to the following:  

```console title="Output"
Launching `main.nf` [condescending_banach] DSL2 - revision: 1bb8f705ad

[93/a56d38] INDEX              | 1 of 1 ✔
[9b/44a943] QUANTIFICATION (1) | 1 of 1 ✔
[gut, /home/ubuntu/hello-nextflow/data/ggal/gut_1.fq, /home/ubuntu/hello-nextflow/data/ggal/gut_2.fq]

```

> Mention .map  

> Mention `row ->` and add graphic

> Mention `file()`

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

!!! question "Poll"

    1. What was the error recieved?  
    2. See the docs on [multiple input channels](https://www.nextflow.io/docs/latest/process.html#multiple-input-channels). Which scope should be amended to fix the issue?
        a. `process { }`
        b. `workflow { }`

    ??? note "Solution" 

        ```console title="Output"
        Process `QUANTIFICATION` declares 2 input channels but 1 were specified
        ```

        ```groovy title="main.nf"
            index_ch = INDEX(params.transcriptome_file)
            quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
        }
        ```

One fixed, your output should look like:  

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

## 2.2.3 Adding a `tag` directive and resuming workflows  

Add `tag` for a more readable execution log. Will also help with profiling
later, when additional samples are added.  

```groovy title="main.nf"
process QUANTIFICATION {

    tag "salmon on ${sample_id}"
    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
    publishDir params.outdir, mode: 'copy'

```

!!! question "Exercise"

    Why would a `tag` directive not be added to `INDEX`?

    ??? note "Solution"

    Only processing once, with a single file.  

Run the workflow with `-resume`:  

```bash
nextflow run main.nf -resume
```

> Mention -resume  

```
Launching `main.nf` [reverent_lavoisier] DSL2 - revision: 10860f201c

executor >  local (2)
[11/010b59] INDEX                          | 1 of 1, cached: 1 ✔
[21/5e9ec8] QUANTIFICATION (salmon on gut) | 1 of 1, cached: 1 ✔
```

> Mention output, cache, resume  

!!! abstract "Summary"

    In this step you have learned:

        1. How to          
        1. How to 
        1. How to 
        1. How to 
        1. How to 

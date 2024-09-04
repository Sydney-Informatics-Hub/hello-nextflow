# 2.2 Samplesheets, operators, and groovy  

!!! note "Learning objectives"

    1. Implement a process with a tuple input
    2. Understand why samplesheets should be used to read in data
    3. Create a channel with operators and Groovy

In this step we will implement `01_fastqc.sh` as a process called `FASTQC` and
use operators and Groovy to create a channel.  

![](img/2.excalidraw.png)

```bash title="01_fastqc.sh"
SAMPLE_ID=gut
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"

mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
fastqc --outdir "results/fastqc_${SAMPLE_ID}_logs" --format fastq ${READS_1} ${READS_2}
```

There's a lot going on in this script, let's break it down.

`SAMPLE_ID=gut` assigns "gut" to the bash variable`SAMPLE_ID`. This is used to:  

- Avoid hardcoding the sample name multiple times in the script  
- Ensure that file pairs of the same sample are processed together  
- Ensure that this script can be run on different sample pairs  

`READS_1` and `READS_2` specify the paths to the `.fq` files.  

Similar to the bash script in the previous step (`00_index.sh`), `mkdir -p`
creates an output folder so that the `fastqc` outputs can be saved here.  

In the `fastqc` command,

- `--outdir` specifies the name of the output directory
- `--format` is a required flag to indicate what format the the reads are in
- `${READS_1}` and `${READS_2}` propagate the paths of the `.fq` files  

## 2.2.1 Building the process  

Start by adding the following `process` scaffold and script definition:  

```groovy title="main.nf"
process FASTQC {
  [ directives ]

  input:
    < process inputs >

  output:
    < process outputs >

  script:
  """
  mkdir -p "fastqc_${sample_id}_logs"
  fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2
  """
}
```

> Note about ${}  

Unlike `salmon` from the previous process, `fastqc` requires that the output
directory be created before running the command, hence the `mkdir ...` line.

!!! question "Exercise"

    Add the output definition

    ??? note "Solution"

        ```groovy title="main.nf"
        output:
        path "fastqc_${sample_id}_logs"
        ```

There are three inputs for this process definition that can be taken from the
script definition you just added:

1. `$sample_id`
2. `$reads_1`
3. `$reads_2`

We will use a tuple as the input qualifier as it's useful to group related
inputs, or, inputs that need to be processed together such as in this case.  

Add the input definition:  

```groovy title="main.nf"
process FASTQC {
  [ directives ]

  input:
  tuple val(sample_id), path(reads_1), path(reads_2)

  output:
  path "fastqc_${sample_id}_logs"

  script:
  """
  fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2
  """
}
```

> Explain tuple a bit more, and the qualifiers i.e. needs brackets

Finish the process by adding the following `container` and `publishDir`
directives. This will be the same in all processes for the workshop, with 
different container values.  

```groovy title="main.nf"
process FASTQC {

    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    publishDir "results", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)
    
    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir -p "fastqc_${sample_id}_logs"
    fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2
    """
}
```

## 2.2.2 Reading files with a samplesheet  

In this step we will prepare the tuple input for the `FASTQC` process using a
samplesheet.

A samplesheet is a delimited text file where each row contains information
or metadata that needs to be processed together.  

https://nextflow-io.github.io/patterns/process-per-csv-record/

Inspect the samplesheet that we will be using:  

```bash
cat data/samplesheet.csv
```

```console title="Output"
sample,fastq_1,fastq_2
gut,data/ggal/gut_1.fq,data/ggal/gut_2.fq
```

The samplesheet has three columns:  

- `sample`: indicates the sample name/prefix
- `fastq_1`, `fastq_2`: contains the relative paths to the reads  

The goal in this step is to read the contents of the samplesheet, and transform
it so it fits the input definition of `FASTQC` we just defined
(`tuple val(sample_id), path(reads_1), path(reads_2)`).

Before that, we need to add an input parameter that points to the samplesheet.  

!!! question "Exercise"

    In your `main.nf` add an input parameter called `reads` and assign the path
    to the samplesheet using `$projectDir`.

    ??? note "Solution"

        ```groovy title="main.nf"
        /*
         * pipeline input parameters
         */
        params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
        params.reads = "$projectDir/data/samplesheet.csv"
        ```

In the next few steps, we will add a mix of Nextflow operators and Groovy
syntax to read in and parse the samplesheet so it is in the correct format
for the process we just added.  

!!! tip

    You do not need to understand what each operator or command does in detail.
    The key takeaway here is to understand that using samplesheets is best
    practice for reading grouped files and metadata into Nextflow, and that
    operators and groovy needs to be chained together to get these in the
    correct format.

Add the following to your workflow scope:  

```groovy title="main.nf"
workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }
        .view()
}
```

> Explain view  

Run the workflow with the `-resume` flag:

```bash
nextflow run main.nf -resume
```

Your output should look something like:  

```console title="Output"
Launching `main.nf` [crazy_einstein] DSL2 - revision: 0ae3776a5e

[de/fef8c4] INDEX [100%] 1 of 1, cached: 1 ✔
[gut, /home/setup2/hello-nextflow/day2/data/ggal/gut_1.fq, /home/setup2/hello-nextflow/day2/data/ggal/gut_2.fq]

```

> Explain -resume, cached  

The chain of commands produces a tuple with three elements that correspond to
the row in the samplesheet. It now fits the requirements of the input
definition of `tuple val(sample_id), path(reads_1), path(reads_2)`, we can
assign the channel to a variable using the
[`set`](https://www.nextflow.io/docs/latest/operator.html#set) operator:

```groovy title="main.nf"
workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
```

We can now call the `FASTQC` process.  

!!! question "Exercise"

    In the `workflow` scope, call the `FASTQC` process with `reads_pairs_ch`
    as the input. Assign it to a variable called `fastqc_ch`

    ??? note "Solution"
    
        ```groovy title="main.nf"
        workflow {
            Channel
                .fromPath(params.reads)
                .splitCsv(header: true)
                .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }
                .set { read_pairs_ch }
        
            index_ch = INDEX(params.transcriptome_file)
            fastqc_ch = FASTQC(read_pairs_ch)
        ```

Run the workflow:  

```bash
nextflow run main.nf -resume
```

Your output should look something like:  

```
Launching `main.nf` [tiny_aryabhata] DSL2 - revision: 9a45f4957b

executor >  local (1)
[de/fef8c4] INDEX      [100%] 1 of 1, cached: 1 ✔
[bb/32a3aa] FASTQC (1) [100%] 1 of 1 ✔
```

If you inspect `results/fastqc_gut_logs` there is an `.html` and `.zip` file
for each of the `.fastq` files.  

??? example "Advanced exercise"  

    Inspect what the `.fromPath()` and `.splitCsv()` commands do by using `.view()`

    ```groovy title="main.nf"
    workflow {
        Channel
            .fromPath(params.reads)
            .view()
    
        index_ch = INDEX(params.transcriptome_file)
    ```
    
    ```console title="Output"
    Launching `main.nf` [hungry_lalande] DSL2 - revision: 587b5b70d1
    
    [de/fef8c4] INDEX [100%] 1 of 1, cached: 1 ✔
    /home/setup2/hello-nextflow/day2/data/samplesheet.csv
    
    ```
    
    ```groovy title="main.nf"
    workflow {
        Channel
            .fromPath(params.reads)
            .splitCsv(header: true)
            .view()
    
        index_ch = INDEX(params.transcriptome_file)
    ```
    
    ```console title="Output"
    Launching `main.nf` [tiny_yonath] DSL2 - revision: 22c2c9d28f
    [de/fef8c4] INDEX | 1 of 1, cached: 1 ✔
    [sample:gut, fastq_1:data/ggal/gut_1.fq, fastq_2:data/ggal/gut_2.fq]
    
    ```

!!! abstract "Summary"

    In this step you have learned:

        1. How to
        1. How to
        1. How to
        1. How to
        1. How to

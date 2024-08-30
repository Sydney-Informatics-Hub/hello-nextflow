# 2.2 Samplesheets, operators, and groovy  

!!! note "Learning objectives"

    1. Implement a process with a tuple input
    2. Understand why samplesheets should be used to read in data
    3. Create a channel with operators and Groovy

In this step we will implement `01_fastqc.sh` as a process called `FASTQC` and
use operators and Groovy to create a channel.  

> Add workflow diagram  

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
  fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2
  """
}
```

> Note about ${}  

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

> Explain tuple a bit more, and the qualifiers  

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
    fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2
    """
}
```

## 2.2.2 Reading files with a samplesheet  

In this step we will prepare the tuple input for the `FASTQC` process using a samplesheet.

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

Add the following lines in the workflow scope of the script:  

```groovy title="main.nf"
workflow {
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)

```

> Add prose about not having to understand in detail but use view to see what
each command does. Add links, things on groovy, operators etc.  

Use some groovy to read and view the samplesheet. Define the channel in the
workflow, and add the
[`view`](https://www.nextflow.io/docs/latest/operator.html#view) operator to see
the contents of the channel:

> Guided example of adding view, then let them do the rest independently  

!!! question "Exercise"

    Add `.view()` after each command (`.fromPath()`, `.splitCsv()`, `.map{}`) and
    run the workflow. Note the output of each step.  

---

## scratch  

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

!!! abstract "Summary"

    In this step you have learned:

        1. How to
        1. How to
        1. How to
        1. How to
        1. How to

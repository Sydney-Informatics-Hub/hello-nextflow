# 2.5 Scaling to multiple samples  

!!! note "Learning objectives"  

    1.

Now that we have a working pipeline on a single-sample, we will update it 
to take multiple samples and introduce Nextflow concepts that help with
understanding and profiling the pipeline.  

## 2.5.1 Labeling tasks with the `tag` directive

The [tag](https://www.nextflow.io/docs/latest/process.html#tag) process
directive allows you to add a custom label, or tag, to each task that gets
executed. It is useful for identifying what is being run when the workflow
is being executed in a bit more detail. It is especially helpful showing you 
what is being run when we run multiple samples, and for profiling later.

Add the following `tag` directives to your existing `FASTQC` and
`QUANTIFICATION` processes. For `FASTQC`:

```groovy title="main.nf" hl_lines="3"
process FASTQC {

    tag "fastqc on ${sample_id}"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    publishDir "results", mode: 'copy'

```

And for `QUANTIFICATION`:  

```groovy title="main.nf" hl_lines="3"
process QUANTIFICATION {

    tag "salmon on ${sample_id}"
    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
    publishDir "results", mode: 'copy'

```

The tags we just added indicates what program is being run (`fastqc` or 
`salmon`), and on which sample (`${sample_id}`) it is being run on. 

We will see this in action in the next section.  

## 2.5.2 Using a samplesheet with multiple samples  

Recall that the samplesheet is used to control which files/data are analysed by
the workflow. Inspect `data/samplesheet_full.csv`.  

```console linenums="1" title="samplesheet_full.csv"
sample,fastq_1,fastq_2
gut,data/ggal/gut_1.fq,data/ggal/gut_2.fq
liver,data/ggal/liver_1.fq,data/ggal/liver_2.fq
lung,data/ggal/lung_1.fq,data/ggal/lung_2.fq
```

Compared to the samplesheet we have been using `data/samplesheet.csv`, this one
contains two additional lines for the `liver` and `lung` paired reads.

Next we will run the workflow with all three samples by inputting 
`data/samplesheet_full.csv` using the double hyphen approach `--` in the run
command.

Run the workflow:  

```bash
nextflow run main.nf -resume --reads data/samplesheet_full.csv
```

Your output should look similar to:  

```console title="Output"
Launching `main.nf` [distraught_bell] DSL2 - revision: dcb06191e7

executor >  local (5)
[de/fef8c4] INDEX                           | 1 of 1, cached: 1 ✔
[4e/b4c797] FASTQC (fastqc on liver)        | 3 of 3, cached: 1 ✔
[36/93c8b4] QUANTIFICATION (salmon on lung) | 3 of 3, cached: 1 ✔
[e7/5d91ea] MULTIQC                         | 1 of 1 ✔

```

There are two new tasks run for `FASTQC` and `QUANTIFICATION`. Our newly added
tags indicate which samples they were run on - either `lung` or `liver` reads!

!!! question "Optional Exercise"

    In the workflow scope, add `.view()` to `reads_in`

    > Update with new workflow definition if keeping this exercise

    ??? note "Solution"

        ```groovy title="main.nf"
        Channel
            .fromPath(params.reads)
            .splitCsv(header: true)
            .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }
            .set { read_pairs_ch }
            read_pairs_ch.view()
        ```

        Run the workflow:  
        
        ```bash
        nextflow run main.nf -resume --reads data/samplesheet_full.csv
        ```
        
        Your output should look something like:  
        
        ```console title="Output"
        
        executor >  local (5)
        [de/fef8c4] INDEX                           | 1 of 1, cached: 1 ✔
        [4e/b4c797] FASTQC (fastqc on liver)        | 3 of 3, cached: 3 ✔
        [36/93c8b4] QUANTIFICATION (salmon on lung) | 3 of 3, cached: 3 ✔
        [e7/5d91ea] MULTIQC                         | 1 of 1 ✔
        [gut, .../data/ggal/gut_1.fq, .../data/ggal/gut_2.fq]
        [liver, .../data/ggal/liver_1.fq, .../data/ggal/liver_2.fq]
        [lung, .../data/ggal/lung_1.fq, .../data/ggal/lung_2.fq]
        ```  

        Key differences to note: 
        
        - Total of three tuples, for each sample  
        - `QUANTIFICATION` and `FASTQC` have 3 processes and 1 cached  
        - Added `results/` outputs for each paired sample  
        - `multiqc_report.html` now has 9 samples  

        Remove `read_pairs_ch.view()` before proceeding.    

> Optional exercise: Move --reads to params.reads in script for downstream nextflow runs?  

## 2.5.3 An introduction to configuration  

> Prose about utilising resources at hand.

> Many ways and things to configure, especially when running on HPC, but beyond
the scope of this workshop. 

We will briefly touch on leveraging multithreading and cpus here.

Some tools like `fastqc` support multithreading. From `fastqc --help`:

```console title="Output"
-t --threads    Specifies the number of files which can be processed    
                simultaneously.
```

Update the `FASTQC` process `script` definition to add this option.  

```groovy title="main.nf"
    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc --outdir "fastqc_${sample_id}_logs" -f fastq $reads_1 $reads_2 -t $task.cpus
    """
```

The [`cpus`](https://www.nextflow.io/docs/latest/process.html#cpus) directive
allows the number of CPUs the process' task should use.  

The `FASTQC` tasks processes paired reads (2 files) per task. Adding
`-t ${task.cpus}` allows them to be processed simultaneously.  

Update your configuration file with the following:  

```groovy title="nextflow.config" hl_lines="1"
process.cpus = 2
docker.enabled = true
```

This will enable 2 cpus to be used per task and will set the `$task.cpus`
in the `FASTQC` process script to use 2 cpus.

## 2.5.4 Profiling  

> Prose about profiling when you add trace, report, timeline to the config.
> This will give insight into task resource requirements and consumption.

To enable the reporting of our workflow execution, add the following to your
`nextflow.config` file:  

```groovy title="nextflow.config" hl_lines="3-5"
process.cpus = 2
docker.enabled = true
trace.enabled = true
report.enabled = true
timeline.enabled = true
```

> Inspect X, Y, Z files by right clicking the file in the VSCode sidebar, and
downloading

!!! abstract "Summary"

    In this step you have learned:

        1. How to
        1. How to
        1. How to

# 2.5 Scaling to multiple samples  

Now that we have a working pipeline on a single-sample, we will update it 
to take multiple samples and introduce Nextflow concepts that help with
understanding and profiling the pipeline.  

## 2.5.1 Adding samples to the samplesheet

Recall that the samplesheet is used to control which files/data are analysed by
the workflow.

!!! question "Exercise"

    Add the `liver` and `lung` samples to the samplesheet by editing
    `data/samplesheet.csv`.  

    ??? tip "Hint"

        Copy the lines where the `gut` samples are defined, and replace the
        sample name with `liver` and `lung`.  

    ??? note "Solution"

        ```bash
        cat data/samplesheet
        ```

        ```console title="Output"
        sample,fastq_1,fastq_2
        gut,data/ggal/gut_1.fq,data/ggal/gut_2.fq
        liver,data/ggal/liver_1.fq,data/ggal/liver_2.fq
        lung,data/ggal/lung_1.fq,data/ggal/lung_2.fq
        ```

Next we will inspect how this changes the workflow.  

!!! question "Exercise"

    In the workflow scope, add `.view()` to see the output of `read_pairs_ch`.

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
nextflow run main.nf -resume
```

Your output should look something like:  

```console title="Output"
Launching `main.nf` [hopeful_boyd] DSL2 - revision: 2d38d1462e

[6e/e2025a] INDEX          [100%] 1 of 1, cached: 1 ✔
[96/81fd93] QUANTIFICATION [100%] 3 of 3, cached: 1 ✔
[5a/f92355] FASTQC         [100%] 3 of 3, cached: 1 ✔
[ab/1eabe7] MULTIQC        [100%] 1 of 1, cached: 1 ✔
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

## 2.5.2 Benchmarking/profiling  

First we need a baseline report of the resource usage per process. The 
`-resume` flag cannot be used here as we need to run the processes again
to record resource usage.  

```bash
nextflow run main.nf -with-report baseline.html
```

- The `-with-report` flag indicates to create an HTML
[execution report](https://www.nextflow.io/docs/latest/tracing.html#execution-report).
- The following argument indicates the output report file name `baseline.html`.  

> Inspect `baseline.html`

Refactor `nextflow.config` and add more cpus per process:   

```groovy linenums="1" title="nextflow.config
process.cpus = 2
docker.enabled = true
```

Run with report:  

```bash
nextflow run main.nf -with-report cpus_2.html
```

> Compare and note the differences (improvements) between reports.  

- run time decreased (8.3s -> 7.7s)  
- allocated cpus of all processes=2
- realtime of all processes decreased in cpus=2
- $cpu in `INDEX` and `FASTQC` increased in cpus=2  

Some programs have options to better utilise resources such as multithreading.

From `fastqc --help`:

```console title="Output"
-t --threads    Specifies the number of files which can be processed    
                simultaneously.
```

`FASTQC` requires multithreading to be explicitly specific in the script.

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

Run again:  

```bash
nextflow run main.nf -with-report fastqc_mt.html
```

> Compare `cpus_2.html` and `fastqc_mt.html` and note the differences between
> %cpu and duration of the FASTQC tasks

- run time decreased  
- %cpu increased  
- duration decreased  

!!! abstract "Summary"

    In this step you have learned:

        1. How to
        1. How to
        1. How to

---

## scratch  

## 2.2.3 Adding a `tag` directive

Add `tag` for a more readable execution log. Will also help with profiling
later, when additional samples are added.

```groovy title="main.nf"
process QUANTIFICATION {

    tag "salmon on ${sample_id}"
    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
    publishDir params.outdir, mode: 'copy'

```

```bash
nextflow run main.nf -resume
```

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
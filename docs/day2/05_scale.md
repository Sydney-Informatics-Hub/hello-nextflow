# 2.5 Scaling to multiple samples  

Next, we will amend the workflow to take in multiple paired-end reads, do some
profiling and optimisation.  

## 2.5.1 Adding samples to the samplesheet

!!! question "Exercise"

    Add the `liver` and `lung` samples to `data/samplesheet.csv`.  

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

!!! question "Exercise"

    Run the workflow. Using `.view()`, what is the new output of `read_pairs_ch`?

    ??? note "Solution"

        ```console title="Output"
        Launching `main.nf` [hopeful_boyd] DSL2 - revision: 2d38d1462e

        [6e/e2025a] INDEX                           [100%] 1 of 1, cached: 1 ✔
        [96/81fd93] QUANTIFICATION (salmon on lung) [100%] 3 of 3, cached: 1 ✔
        [5a/f92355] FASTQC (fastqc on liver)        [100%] 3 of 3, cached: 1 ✔
        [ab/1eabe7] MULTIQC                         [100%] 1 of 1, cached: 1 ✔
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

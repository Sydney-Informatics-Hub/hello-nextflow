# 2.5 Scaling to multiple samples  

!!! note "Learning objectives"  

    1.

Now that we have a working pipeline on a single-sample, we will update it 
to take multiple samples and introduce Nextflow concepts that help with
understanding and profiling the pipeline.  

## 2.5.1 Using a samplesheet with multiple samples  

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

- Two new tasks run for FASTQC and QUANTIFICATION 
- Tags show different sample (non-gut)

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

> Remove `read_pairs_ch.view()` before proceeding.    

> Move --reads to params.reads in script for downstream nextflow runs?  

## 2.5.2 Introduction to configuration  

Prose about utilising resources at hand.

Many ways and things to configure, especially when running on HPC, but beyond
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

Update your configuration file:  

```groovy linenums="1" title="nextflow.config
process.cpus = 2
docker.enabled = true
```

> Show what happens when you add trace, report, timeline to the config.
This will give insight into task resource requirements and consumption.

```groovy linenums="1" title="nextflow.config
process.cpus = 2
docker.enabled = true
trace.enabled = true
report.enabled = true
timeline.enabled = true
```

!!! abstract "Summary"

    In this step you have learned:

        1. How to
        1. How to
        1. How to

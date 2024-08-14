# 2.4 Scaling to multiple samples  

Next, we will amend the workflow to take in multiple paired-end reads, do some
profiling and optimisation.  

## 2.4.1 Using glob to take multiple samples  

> What parameter should be changed to take in additional paired-reads?  

> What should the updated parameter be amended to? Tip: see [docs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)  

First, update `params.reads`:  
```
params.reads = "$projectDir/data/ggal/*_{1,2}.fq"
```

> View the output of the channel

```
workflow {
    ...
   
    read_pairs_ch.view()
}
```

Run with `-resume`.  

```
$ nextflow run main.nf -resume  

...
executor >  local (5)
[98/bc54e7] INDEX                            | 1 of 1, cached: 1 ✔
[13/b79c22] QUANTIFICATION (salmon on lung)  | 3 of 3, cached: 1 ✔
[42/3c2765] FASTQC (FASTQC on liver)         | 3 of 3, cached: 1 ✔
[9a/758ecc] MULTIQC                          | 1 of 1 ✔
[lung, [/home/temp/GitHub/hello-nextflow/data/ggal/lung_1.fq, /home/temp/GitHub/hello-nextflow/data/ggal/lung_2.fq]]
[gut, [/home/temp/GitHub/hello-nextflow/data/ggal/gut_1.fq, /home/temp/GitHub/hello-nextflow/data/ggal/gut_2.fq]]
[liver, [/home/temp/GitHub/hello-nextflow/data/ggal/liver_1.fq, /home/temp/GitHub/hello-nextflow/data/ggal/liver_2.fq]]

```

Key differences to note: 
- Total of three tuples, for each sample  
- `QUANTIFICATION` and `FASTQC` have 3 processes and 1 cached  
- Added `results/` outputs for each paired sample  
- `multiqc_report.html` now has 9 samples  

> Remove .view() channelfor next steps. 

## 2.4.2 Profiling/scaling/optimisation  

First we need a baseline report of the resource usage per process. The 
`-resume` flag cannot be used here as we need to run the processes again
to record resource usage.  

```
$ nextflow run main.nf -with-report baseline.html
```

- The `-with-report` flag indicates to create an HTML
[execution report](https://www.nextflow.io/docs/latest/tracing.html#execution-report).
- The following argument indicates the output report file name `baseline.html`.  

> Inspect `baseline.html`

Refactor `nextflow.config` and add more cpus per process:   

```
process {
    cpus = 2
    container = 'nextflow/rnaseq-nf'
}

docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}
```

Run with report:  

```
nextflow run main.nf -with-report cpus_2.html
```

> Compare and note the differences (improvements) between reports.  

- run time decreased (8.3s -> 7.7s)  
- allocated cpus of all processes=2
- realtime of all processes decreased in cpus=2
- $cpu in `INDEX` and `FASTQC` increased in cpus=2  

Some programs have options to better utilise resources such as multithreading.

Add multithreading for `FASTQC`:  

```
script:
"""
mkdir fastqc_${sample_id}_logs
fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads} -t ${task.cpus}
"""
```

Run again:  

```
nextflow run main.nf -with-report fastqc_mt.html
```
From `fastqc --help`:

```
-t --threads    Specifies the number of files which can be processed    
                simultaneously.
...
```

`FASTQC` processes read pairs (2 files) per task. The `-t` flag allows them to
be processed at the same time using `${task.cpus}` (2 cpus defined in .config).  

> Compare `cpus_2.html` and `fastqc_mt.html` and note the differences between
> %cpu and duration of the FASTQC tasks

- run time decreased  
- %cpu increased  
- duration decreased  

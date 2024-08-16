# 2.4 MultiQC report

This step collects the outputs from the `QUANTIFICATION` and `FASTQC` processes
to create a final report using `multiqc`.  

> collect, reduce graphics/description  

## 2.4.1 Combining channels with operators  

This is a snippet of the help page for `multiqc`:  

```console
Usage: multiqc [OPTIONS] <analysis directory>

  MultiQC aggregates results from bioinformatics analyses across many
  samples into a single report.

  It searches a given directory for analysis logs and compiles a HTML
  report. It's a general use tool, perfect for summarising the output from
  numerous bioinformatics tools.

  To run, supply with one or more directory to scan for analysis results. To
  run here, use 'multiqc .'
```

Importantly, it indicates that the inputs are "one or more directories".  

In our example, the inputs will be `fastqc_gut_logs` and `gut`.  

!!! question "Exercise"

    Which channels output these?

    ??? note "Solution"

        `fastqc_ch` and `quant_ch`.  

The next few steps will involve chaining together nextflow operators to prepare
the inputs for the correct format for the `MULTIQC` process.  

In the workflow scope, use the 
[`mix`](https://www.nextflow.io/docs/latest/operator.html#mix) operator to
emit the contents of `fastqc_ch` and `quant_ch` and view it:  

```groovy title="main.nf"
    fastqc_ch = FASTQC(read_pairs_ch)

    quant_ch
        .mix(fastqc_ch)
        .view()
}
```

Run the workflow:  

```bash
nextflow run main.nf -resume  
```

The output should look something like:  

```console title="Output"
Launching `main.nf` [hopeful_perlman] DSL2 - revision: c0603d2f16

[aa/3b8821] INDEX                          | 1 of 1, cached: 1 ✔
[c2/baa069] QUANTIFICATION (salmon on gut) | 1 of 1, cached: 1 ✔
[ad/e49b20] FASTQC (fastqc on gut)         | 1 of 1, cached: 1 ✔
/home/ubuntu/hello-nextflow/work/ad/e49b204c50c8f819dbddc4526010cf/fastqc_gut_logs
/home/ubuntu/hello-nextflow/work/c2/baa0691de7afccb3ed6d285947c8bf/gut

```  

The outputs have been emitted one after the other, meaning that it will be
processed separately. 

Add the [`collect`](https://www.nextflow.io/docs/latest/operator.html#collect)
operator to ensure all samples are processed together in the same
process and view the output:  

```groovy title="main.nf"
    fastqc_ch = FASTQC(read_pairs_ch)

    quant_ch
        .mix(fastqc_ch)
        .collect()
        .view()
}
```

Run the workflow:  

```bash
nextflow run main.nf -resume  
```

The channel now outputs a single tuple with the two directories.  

```console title="Output"
Launching `main.nf` [cranky_booth] DSL2 - revision: 568fddbfec

[aa/3b8821] INDEX                          | 1 of 1, cached: 1 ✔
[c2/baa069] QUANTIFICATION (salmon on gut) | 1 of 1, cached: 1 ✔
[ad/e49b20] FASTQC (FASTQC on gut)         | 1 of 1, cached: 1 ✔
[/home/temp/GitHub/hello-nextflow/work/ad/e49b204c50c8f819dbddc4526010cf/fastqc_gut_logs, /h
ome/temp/GitHub/hello-nextflow/work/c2/baa0691de7afccb3ed6d285947c8bf/gut]

```

!!! question "Exercise"

    Remove `.view()` and `set` the output as `all_qc_ch` 

    ??? note "Solution"
        
        ```groovy title="main.nf"
            quant_ch
                .mix(fastqc_ch)
                .collect()
                .set { all_qc_ch }
        }
        ```

Add the following process definition:  

```groovy title="main.nf"
process MULTIQC {
   
   container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
   publishDir params.outdir, mode: 'copy'

   input:
   path '*'

   output:
   path 'multiqc_report.html'

   script:
   """
   multiqc .
   """
}
```

!!! question "Exercise"

    Add the workflow scope for `MULTIQC`.

    ??? note "Solution"

```groovy title="main.nf"
    quant_ch
        .mix(fastqc_ch)
        .collect()
        .set { all_qc_ch }

    MULTIQC(all_qc_ch)
}
```

Run the workflow:  

```bash
nextflow run main.nf -resume  
```

```console title="Output"
Launching `main.nf` [hopeful_swanson] DSL2 - revision: a4304bbe73

[aa/3b8821] INDEX                          [100%] 1 of 1, cached: 1 ✔
[c2/baa069] QUANTIFICATION (salmon on gut) [100%] 1 of 1, cached: 1 ✔
[ad/e49b20] FASTQC (FASTQC on gut)         [100%] 1 of 1, cached: 1 ✔
[a3/1f885c] MULTIQC                        [100%] 1 of 1 ✔

```

> Inspect `results/multiqc_report.html`, maybe Poll something in the file  

You have a working pipeline for a single paired-end sample!

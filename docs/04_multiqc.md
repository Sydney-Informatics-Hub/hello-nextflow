## 2.4 MultiQC report

This step collects the outputs from the `QUANTIFICATION` and `FASTQC` processes
to create a final report using `multiqc`.  

### 2.4.1 Implementing the process  

This is a snippet of the help page for `multiqc`.  

```
Usage: multiqc [OPTIONS] <analysis directory>

  MultiQC aggregates results from bioinformatics analyses across many
  samples into a single report.

  It searches a given directory for analysis logs and compiles a HTML
  report. It's a general use tool, perfect for summarising the output from
  numerous bioinformatics tools.

  To run, supply with one or more directory to scan for analysis results. To
  run here, use 'multiqc .'
```

```groovy title="main.nf"
process MULTIQC {

    script:  
    """
    multiqc .
    """
}
```

The inputs will be `fastqc_gut_logs` and `gut`.  

!!! question "Exercise"

    Which channels output these?

    ??? note "Solution"

        `fastqc_ch` and `quant_ch`.  

```
process MULTIQC {
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

### 2.4.2 Creating a channel with operators  

Before adding the `MULTIQC` process in the workflow, we first need to figure
out how to pass the QC directories as input.  

Explanation of Gather - reduce - set/emit.  

Use the `mix` and `collect` operators to gather the outputs from
`QUANTIFICATION` and `FASTQC`  

```
workflow {
    ...

    quant_ch
        .mix(fastqc_ch)
        .view()
}
```

```
$ nextflow run main.nf -resume  

Launching `06.nf` [hopeful_perlman] DSL2 - revision: c0603d2f16

[aa/3b8821] INDEX                          | 1 of 1, cached: 1 ✔
[c2/baa069] QUANTIFICATION (salmon on gut) | 1 of 1, cached: 1 ✔
[ad/e49b20] FASTQC (FASTQC on gut)         | 1 of 1, cached: 1 ✔
/home/temp/GitHub/hello-nextflow/work/ad/e49b204c50c8f819dbddc4526010cf/fastqc_gut_logs
/home/temp/GitHub/hello-nextflow/work/c2/baa0691de7afccb3ed6d285947c8bf/gut

```  

Then use `collect` to ensure all samples are processed together in the same
process and view the output:  
```
workflow {
    ...

    quant_ch
        .mix(fastqc_ch)
        .collect()
        .view()
}
```

Run.  

```
Launching `main.nf` [cranky_booth] DSL2 - revision: 568fddbfec

R N A S E Q - N F   P I P E L I N E
===================================
transcriptome: /home/temp/GitHub/hello-nextflow/data/ggal/transcriptome.fa
reads        : /home/temp/GitHub/hello-nextflow/data/ggal/gut_{1,2}.fq
outdir       : results

[aa/3b8821] INDEX                          | 1 of 1, cached: 1 ✔
[c2/baa069] QUANTIFICATION (salmon on gut) | 1 of 1, cached: 1 ✔
[ad/e49b20] FASTQC (FASTQC on gut)         | 1 of 1, cached: 1 ✔
[/home/temp/GitHub/hello-nextflow/work/ad/e49b204c50c8f819dbddc4526010cf/fastqc_gut_logs, /h
ome/temp/GitHub/hello-nextflow/work/c2/baa0691de7afccb3ed6d285947c8bf/gut]

```

Important to note that the channel now emits a tuple.  

Now set the output to `allqc_ch` and input to `MULTIQC`.  

```
workflow {
    ...

    quant_ch
        .mix(fastqc_ch)
        .collect()
        .set { allqc_ch }

    MULTIQC(allqc_ch)
}
```

Run with the `-resume` flag. You should have a `results/multiqc_report.html`.
View the file in a browser - should have 3 samples, a section for `salmon` and
`fastqc` with some QC stats etc.  

You have a working pipeline for a single paired-end sample!

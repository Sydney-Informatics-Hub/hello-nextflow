# 2.2 Perform expression quantification

## 2.2.1 Collect read files by pairs  

There are numerous channel factories that can be used to create channels. In
this step, you will use the
[fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)
channel factory to create a channel of read pairs.  

The `fromFilePairs` channel factory takes a glob pattern as input and returns
a channel of tuples.

Add a reads parameter that takes as input the paired gut reads:  

```groovy title="main.nf"
/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.outdir = "results"
```

!!! question "Exercise"

    Add `params.reads` to `log.info`  

    ??? Solution

        ```groovy title="main.nf"
        log.info """\
            R N A S E Q - N F   P I P E L I N E
            ===================================
            transcriptome: ${params.transcriptome_file}
            reads        : ${params.reads}
            outdir       : ${params.outdir}
            """
            .stripIndent()
        ```

Need minimal explanation of "factories" and "tuples" - use graphic here.

i.e. groups by the shared pattern/prefix  

Add links to docs/examples for tuples.  

Define the channel in the workflow with the `.view()` command to see the
contents of the channel:  

```groovy title="main.nf"
workflow {
    Channel
        .fromFilePairs(params.reads)
        .view()

    index_ch = INDEX(params.transcriptome_file)
```

Run the workflow.

```bash
nextflow run main.nf  
```

Your output will display a tuple consisting of two elements.  

```console title="Output"
Launching `02.nf` [drunk_descartes] DSL2 - revision: 4698cdd674

R N A S E Q - N F   P I P E L I N E
===================================
transcriptome: /home/fredjaya/GitHub/hello-nextflow/data/ggal/transcriptome.fa
reads        : /home/fredjaya/GitHub/hello-nextflow/data/ggal/gut_{1,2}.fq
outdir       : results

[-        ] INDEX | 0 of 1
executor >  local (1)
executor >  local (1)
[7b/e0040f] INDEX | 1 of 1 ✔
[gut, [/home/fredjaya/GitHub/hello-nextflow/data/ggal/gut_1.fq, /home/fredjaya/GitHub/hello-nextflow/data/ggal/gut_2.fq]]

```

The first element of the tuple (`gut`) is the read pair prefix, and the second
is a list representing the files.  

We will explore how glob patterns can be used in later steps.  

The [`set`](https://www.nextflow.io/docs/latest/operator.html#set) operator can
also be used to define a new channel variable in place
of an `=` assignment.  

Remove `.view()` and assign the channel output with `set`:  

```groovy title="main.nf"
workflow {
    Channel
        .fromFilePairs(params.reads)
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
```


## 2.2.2 Implementing the process  

Add the following process definition. 
```
process QUANTIFICATION {

	input:
	path salmon_index
	tuple val(sample_id), path(reads)

	output:
	path "$sample_id"

	script:
	"""
	salmon quant --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
	"""
```

> Add the `publishDir` directive to the process.  
```
process QUANTIFICATION {
	publishDir params.outdir, mode: 'copy'

	input:
	...
```

> Update the workflow, assigning the output of `QUANTIFICATION` to `quant_ch` 

https://www.nextflow.io/docs/latest/process.html#multiple-input-channels
```
workflow {
    Channel
		.fromFilePairs(params.reads)
		.set { read_pairs_ch }

	index_ch = INDEX(params.transcriptome_file)
	quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)

```

Add `tag` for a more readable execution log. Will also help with profiling
later, when additional samples are added.  

> Question: Why would a `tag` directive not be added to `INDEX`?  

Run the workflow.
```
Launching `main.nf` [reverent_lavoisier] DSL2 - revision: 10860f201c

R N A S E Q - N F   P I P E L I N E
===================================
transcriptome: /home/fredjaya/GitHub/hello-nextflow/data/ggal/transcriptome.fa
reads        : /home/fredjaya/GitHub/hello-nextflow/data/ggal/gut_{1,2}.fq
outdir       : results

executor >  local (2)
[11/010b59] INDEX                          | 1 of 1 ✔
[21/5e9ec8] QUANTIFICATION (salmon on gut) | 1 of 1 ✔
```

> Question: Check `results/`, what is the name of the output directory?  

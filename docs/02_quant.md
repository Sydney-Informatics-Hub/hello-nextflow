# 2.2 Perform expression quantification

## 2.2.1 Collect read files by pairs  

Add `params.reads` to input parameters.

Mention channels, and `fromFilePairs`.  

The `fromFilePairs` channel factory takes a glob pattern as input and returns
a channel of tuples.

> Add `params.reads` to `log.info`  

```
/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.outdir = "results"

log.info """\
    R N A S E Q - N F   P I P E L I N E
	===================================
	transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
	outdir       : ${params.outdir}
	"""
	.stripIndent(true)

```

Need minimal explanation of "factories" and "tuples" - use graphic here.

i.e. groups by the shared pattern

Define the channel in the workflow. We will also add `.view()` to see the 
contents:
```
workflow {
	Channel
		.fromFilePairs(params.reads)
		.set { read_pairs_ch }
	read_pairs_ch.view()

	index_ch = INDEX(params.transcriptome_file)
```

Mention  `.set` and `.view`.  

Remove `read_pairs_ch.view()` and run.

```
$ nextflow run main.nf  

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

Add links to docs/examples for tuples.  

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

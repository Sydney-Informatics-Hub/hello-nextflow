# Day 2: Simple RNA-Seq workflow  

Introduction...  

## 2.1 Indexing a transcriptome file

The first step in the RNA-Seq workflow is to index the transcriptome
using `salmon`. The minimum command required is:  

```bash
salmon index -t [transcriptome_file] -i salmon_index
```

Where:  
- `-t` is the flag for the transcriptome file (`data/ggal/transcriptome.fa`).  
- `-i` is the name of the directory where `salmon` index files will be output.  

Note that `-i` outputs the directory and does not need to be piped in the
process script definition.

1. Define `params.transcriptome_file`  

`projectDir` indicates where the main script is located.  

```
/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
```

2. Define a process called INDEX (input, output, script)
```
/*
 * define the `INDEX` process that create a binary index
 * given the transcriptome file
 */
process INDEX {
	input:
	path transcriptome

	output:
	path 'salmon_index'

	script:
	"""
	salmon index -t $transcriptome -i salmon_index
	"""
```

3. Add the workflow scope  
```
workflow {
    index_ch = INDEX(params.transcriptome_file)
}
```

Run nextflow
```
$ nextflow run main.nf

N E X T F L O W   ~  version 24.04.3

Launching `main.nf` [awesome_kimura] DSL2 - revision: 2d008f1c4f

executor >  local (1)
[9a/4a1dc7] INDEX | 1 of 1 ✔
```

Inspect files, `salmon_index` contains a bunch of files etc..

Some QoL changes.  

In Day 1, we had to inspect the `workDir`, add `publishDir` to `INDEX` to output in a human-readable and central location. Add `params.outdir` and `publishDir`:
```
/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.outdir = "results"
```

```
process INDEX {
	publishDir params.outdir, mode: 'copy'

	input:
```

Explain `mode: 'copy'` and link to docs.  

Explain groovy.  

Add a `log.info` to track information about the workflow, mainly parameters.  
```
log.info """\
    R N A S E Q - N F   P I P E L I N E
	===================================
	transcriptome: ${params.transcriptome_file}
	outdir       : ${params.outdir}
	"""
	.stripIndent(true)
```

Explain `.stripIndent()`.  

```
$ nextflow run main.nf

N E X T F L O W   ~  version 24.04.3

Launching `main.nf` [goofy_ekeblad] DSL2 - revision: 8420244f03

R N A S E Q - N F   P I P E L I N E
===================================
transcriptome: /home/fredjaya/GitHub/hello-nextflow/data/ggal/transcriptome.fa
outdir       : results

executor >  local (1)
[c7/83b9f9] INDEX | 1 of 1 ✔

```

`log.info` prints the parameters.  

Now all the index files are in `/results` folder, we will add all outputs for
the remaining processes here.

## 2.2 Perform expression quantification

### 2.2.1 Collect read files by pairs  

Add `params.reads` to input parameters and `log.info`:  

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

Explanation about channels, and `fromFilePairs`.  

The `fromFilePairs` channel factory takes a glob pattern as input and returns
a channel of tuples.

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

Explanation about `.set` and `.view`.  

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


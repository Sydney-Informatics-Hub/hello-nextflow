# Day 2: Simple RNA-Seq workflow  

Introduction...  

First make a working workflow for a single paired-end sample.

Then scale for multiple samples and optimise.

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

### 2.1.1 Implementation  

1. Define `params.transcriptome_file`  

`projectDir` indicates where the main script is located.  

```
/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
```

> Given the process script definition, complete the following:
> 	- `input`
>	- `output`

Define a process called INDEX (input, output, script)
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

> Complete the workflow scope so the output for `INDEX` is assigned to a
> channel named `index_ch`

Add the workflow scope  
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

### 2.1.2 Nextflow "housekeeping"  

Now that we have a single process working, we will add a few things to tidy the
workflow and the ouputs.  

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

Mention `mode: 'copy'` and link to docs.  

Mention groovy.  

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

**Achievements**  

## 2.2 Perform expression quantification

### 2.2.1 Collect read files by pairs  

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

### 2.2.2 Implementing the process  

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

## 2.3 Quality control

Next, you will implement a `FASTQC` quality control step for your input reads.
The inputs are the same as the read pairs used in the `QUANTIFICATION` step.

### 2.3.1 Implementation  

Add the following to the workflow:  
```
workflow {
	...
	fastqc_ch = FASTQC(read_pairs_ch)
}
```

Add the following `FASTQC` process.

> Complete the process by adding the `input`, `tag`, and `publishDir`
> directives

```
process FASTQC {
    tag "FASTQC on $sample_id" // leave empty
	publishDir $params.outdir, mode: 'copy'
	input:
	tuple val(sample_id), path(reads) // leave empty

	output:
	path "fastqc_${sample_id}_logs"

	script:
	"""
	mkdir fastqc_${sample_id}_logs
	fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
	"""
}
```

Run the workflow. 

> What is output?  

```
...

Comand error:
 .command.sh: line 3: fastqc: command not found

...
```

### 2.3.2 Using containers  

This execution will fail because `fastqc` is not installed in your environment.
A docker (**to be replaced with Singularity/Apptainer**) container image with
the remaining software required will be used. 

Create a file called `nextflow.config` in the same directory as `main.nf`. Add
the following and re-run the workflow:
```
process.container = 'nextflow/rnaseq-nf'
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g)'
```

This time the execution will work because `docker.enabled = true` tells
nextflow to use the docker container `nextflow/rnaseq-nf`.

Check that the output files have been output correctly in `/results`.  

## 2.4 MultiQC report

This step colelcts the outputs from the `QUANTIFICATION` and `FASTQC` processes
to create a final report using `multiqc`.  

### 2.4.1 Creating a channel with operators  

Gather - reduce - set/emit.  

Use the `mix` and `collect` operators to gather the outputs from
`QUANTIFICATION` and `FASTQC`  

---  

**To-dos**  
- [ ] Introduce `--resume` early on somewhere  
- [ ] Make a separate input channel for `QUANTIFICATION`  
- [ ] Replace docker with singularity/apptainer  
- [ ] Stretch tasks per section  

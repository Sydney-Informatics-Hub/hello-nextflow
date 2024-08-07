# 2.1 Indexing a transcriptome file

The first step in the RNA-Seq workflow is to index the transcriptome
using `salmon`. The minimum command required is:  

```bash
salmon index -t [transcriptome_file] -i salmon_index
```

Where:  
- `-t` is the flag for the transcriptome file (`data/ggal/transcriptome.fa`).  
- `-i` is the name of the directory where `salmon` index files will be output.  

Note that `-i` outputs the directory `salmon_index` and does not need to be piped in the
process script definition.

## 2.1.1 Add the `INDEX` process

In the empty `main.nf` script, define the  `params.transcriptome_file`:  

```groovy title="main.nf"
/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
```

`$projectDir` indicates where the `main.nf` script is located.

> Using the information provided in the previous section, complete the process definition that takes in a single input for the transcriptome. Ensure that the process includes definitions for the `input`, `output`, and `script`.  

<details><summary>Show code:</summary>
<br>
```groovy
/*
 * define the `INDEX` process that creates a binary index
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
</details>

---

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

## 2.1.2 Nextflow "housekeeping"  

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

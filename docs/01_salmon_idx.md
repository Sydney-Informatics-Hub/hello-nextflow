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

In the empty `main.nf` script, define the `params.transcriptome_file`:  

```groovy linenums="1" title="main.nf"
/*
 * pipeline input parameters
 */
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
```

Recall that parameters are inputs and options that can be modified when the
workflow is executed.  

`$projectDir` indicates the directory the `main.nf` script is located.

!!! question "Exercise"

    Using the information provided in the previous section, complete the process
    definition that takes in a single input for the transcriptome. Ensure that
    the process includes definitions for the `input`, `output`, and `script`.  

    ??? Solution

        ```groovy linenums="6" title="main.nf"
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
        }
        ```

!!! question "Exercise"

    Complete the workflow scope so the output for `INDEX` is assigned to a
    channel named `index_ch`

    ??? Solution

        ```groovy title="main.nf"
        workflow {
            index_ch = INDEX(params.transcriptome_file)
        }
        ```

Run the workflow using the following command:  

```bash
nextflow run main.nf
```  

Your output will look something like this:  

```console title="Output"
N E X T F L O W   ~  version 24.04.3

Launching `main.nf` [awesome_kimura] DSL2 - revision: 2d008f1c4f

executor >  local (1)
[9a/4a1dc7] INDEX | 1 of 1 ✔
```

Recall that the specifics of the output are randomly generated (i.e.
`[awesome_kimura]` and `[9a/4a1dc7]` in this example).

In this example, the output files for the `INDEX` process is output in
`work/9a/4a1dc7...`. As more processes are added, it becomes cumbersome to
navigate and inspect each directory for the output files.  

## 2.1.2 Save files to an output directory with `publishDir`  

Next, you will create a centralised directory where all output files will be
copied to with the `publishDir` directive.  

!!! question "Exercise"

    Add a second parameter named `outdir`, and assign the string `"results"`.

    ??? Solution

        ```groovy linenums="1" title="main.nf"
        /*
         * pipeline input parameters
         */
        params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
        params.outdir = "results"
        ```

Add the `publishDir` directive to `INDEX`:  

```groovy title="main.nf"
process INDEX {
    publishDir params.outdir, mode: 'copy'

```

`mode: 'copy'` copies the output files into the publish directory (`results/`).
More information and other modes can be found on
[publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir).

## 2.1.3 Logging information with `log.info`  

It can be useful to print the parameters when workflows are executed.

The `log.info` command can be used to print multiline information using
[groovy](https://www.tutorialspoint.com/groovy/groovy_basic_syntax.htm)'s
logger functionality.

Add the `log.info` command to print the parameters:  
```groovy title="main.nf"
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    outdir       : ${params.outdir}
    """
    .stripIndent()
```

`.stripIndent()` is a groovy method that removes the indentation at the
beginning of the lines. We will use groovy more in later exercises.

Run the workflow:  

```bash
nextflow run main.nf
```

```console title="Output"
N E X T F L O W   ~  version 24.04.3

Launching `main.nf` [goofy_ekeblad] DSL2 - revision: 8420244f03

R N A S E Q - N F   P I P E L I N E
===================================
transcriptome: /home/fredjaya/GitHub/hello-nextflow/data/ggal/transcriptome.fa
outdir       : results

executor >  local (1)
[c7/83b9f9] INDEX | 1 of 1 ✔

```

Your output now displays the parameter, and the output index files are in
`results/salmon_index`.  

All remaining output files generated with the workflow will now be copied from
the `workDir` to here.

!!! abstract "Summary"

    In this step you have learned:  

    1. How to implement a simple process with input data  
    2. How to define parameters in your workflow scripts  
    3. How to improve workflows with Nextflow directives and groovy  

# 2.1 Indexing a transcriptome file

The first step in the RNA-Seq workflow is to index the transcriptome
using `salmon`. The minimum command required is:  

```bash
salmon index -t [transcriptome_file] -i salmon_index
```

Where:  

- `-t` is the flag for the transcriptome file (`data/ggal/transcriptome.fa`).  
- `-i` is the name of the directory where `salmon` index files will be output.  

Note that `-i` outputs the directory `salmon_index` and does not need to be 
piped in the process script definition.

(Usually should run manually on CLI to test I/O but introduces docker
complexity)  

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

    Using the information provided in the previous section, define a process
    called `INDEX`. The process takes in a single input for the transcriptome
    file. Ensure that the process includes definitions for the `input`,
    `output`, and `script`.  

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

Your output will return an error that looks like:  

```console title="Output"
N E X T F L O W   ~  version 24.04.4

Launching `main.nf` [awesome_kimura] DSL2 - revision: 2d008f1c4f

executor >  local (1)
[9a/4a1dc7] INDEX | 1 of 1 ✔
ERROR ~ Error executing process > 'INDEX'

Caused by:                                    
  Process `INDEX` terminated with an error exit status (127)


Command executed:

  salmon index -t transcriptome.fa -i salmon_index

Command exit status:
  127

Command output:
  (empty)

Command error:                                
  .command.sh: line 2: salmon: command not found

Work dir:                                     
  /home/ubuntu/hello-nextflow/work/f9/75c584fde0f36c40431b4547318bb7

Tip: you can replicate the issue by changing to the process work dir and entering the comma
nd `bash .command.run` 

 -- Check '.nextflow.log' file for details

```

The error above is caused because `salmon` is not installed in your environment

## 2.1.2 Using docker containers  

Nextflow has support for managing the execution of processes in Docker
containers. This is useful when you need to execute a process that requires a
specific software version.  

Add the following line to `INDEX`:  

```groovy title="main.nf"
process INDEX {

    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"

    input:
```

Usually, containers need to be downloaded using a command such as
`docker pull [image]'. All containers have been saved in saved in the
environment beforehand.  

Run the workflow:  

```bash
nextflow run main.nf
```  

You will see that a similar output is returned. Execution with containers must
be explicitly specified.  

Create a file called `nextflow.config` in the same directory as `main.nf`.  

Add the following line:  

```groovy linenums="1" title="nextflow.config"
docker.enabled = true
```

And run the workflow:  

```bash
nextflow run main.nf
```  

The output should like something like:  

```console title="Output"
N E X T F L O W   ~  version 24.04.4

Launching `01_b.nf` [loquacious_elion] DSL2 - revision: 3e244cb317

executor >  local (1)
[26/c410b1] INDEX | 1 of 1 ✔

```

Recall that the specifics of the output are randomly generated (i.e.
`[loquacious_elion]` and `[26/c410b1]` in this example).

In this example, the output files for the `INDEX` process is output in
`work/26/c410b1...`. As more processes are added, it becomes cumbersome to
navigate and inspect each directory for the output files.  

## 2.1.3 Save files to an output directory with `publishDir`  

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
    
    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
    publishDir params.outdir, mode: 'copy'

    input:
```

`mode: 'copy'` copies the output files into the publish directory (`results/`).
More information and other modes can be found on
[publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir).

Run the workflow:  

```bash
nextflow run main.nf
```  

The files are now output to `results/salmon_index`.  

`container` and `publishDir` will be included in each following process. The
remaining output files generated with the workflow will now be copied from the
`workDir` to here.

!!! abstract "Summary"

    In this step you have learned:  

        1. How to implement a simple process with input data  
        2. How to define parameters in your workflow scripts  
        3. How to use a docker container for a process  
        4. How to output files in a dedicated `publishDir`  

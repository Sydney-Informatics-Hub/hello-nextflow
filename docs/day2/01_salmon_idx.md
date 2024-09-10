# 2.1 Implementing a simple process with Docker  

!!! info "Learning objectives"

    1. Implement a Nextflow process that takes a single file as input.  
    2. Understand the importance of containers in ensuring consistent and
    reproducible execution across processes.
    3. Store output files with the `publishDir` directive.  

In this step we will implement `00_index.sh` as a Nextflow process called
`INDEX`. 

![](img/1.excalidraw.png)

```bash title="00_index.sh"
mkdir "results"
salmon index \
    --transcripts data/ggal/transcriptome.fa \
    --index results/salmon_index
```

- The script first creates a `results/` folder then runs `salmon index`.  
- The `--transcripts` flag indicates that the path to the input transcriptome
file is (`data/ggal/transcriptome.fa`).  
- `--index results/salmon_index` tells `salmon` to save the output index files in a directory called `salmon_index`, within the newly created `results` directory.  

> Note about hardcoded arguments and output directory, and how this is
addressed in Nextflow  

This is enough information to construct the process.  

## 2.1.1 Adding the `INDEX` process

In the empty `main.nf` script, add the following `process` scaffold with the
script definition:  

```groovy title="main.nf"
process INDEX {
  [ directives ]

  input:
    < process inputs >

  output:
    < process outputs >

  script:
  """
  salmon index --transcripts $transcriptome --index salmon_index
  """
}
```

Next, we will edit the `input` and `output` definitions to match the specific
data and results for this process. In the `00_index.sh` script, the relevant
information is:  

* The fasta (`.fa`) file defined by the variable `$transcriptome` and provided
to the `--transcripts` flag  
* The index output directory `salmon_index/` provided to the `--index` flag  

!!! info "Defining inputs and outputs"

    Remember, input and output definitions require a qualifier and name. For example:  
    ```groovy
    input:
    <input qualifier> <input name>

    output:
    <output qualifier> <output name>
    ```

    The qualifier defines the type of data, and the names are treated like variables.  

```groovy title="main.nf"
process INDEX {
  [ directives ]

  input:
  path transcriptome

  output:
  path 'salmon_index'

  script:
  """
  salmon index --transcripts $transcriptome --index salmon_index
  """
}
```

Note that the input `path transcriptome` refers to a variable, meaning the
actual file or directory provided as input can be changed depending on the data
you provide it. The output `path 'salmon_index'` is fixed, meaning it will
always create an output folder called `salmon_index`, no matter what the input
is.  

This is how Nextflow can handle different inputs while always producing the
same output name.  

More information on using input and output blocks can be found in the process
[inputs](https://www.nextflow.io/docs/latest/process.html#inputs) and
[outputs](https://www.nextflow.io/docs/latest/process.html#outputs) Nextflow
documentation, respectively.  

> Add your own comment 

## 2.1.2 Save files to an output directory with `publishDir`  

Next we will implement the Nextflow equivalent of saving the output files into a
`results/` directory.  

Replace `[ directives ]` in your `main.nf` script with the `publishDir` 
directive, specifying the directory name as `"results"` and the mode as
`'copy'`. Your `main.nf` should look like this: 

```groovy title="main.nf"
process INDEX {

  publishDir "results", mode: 'copy'

  input:
  path transcriptome

  output:
  path 'salmon_index'

  script:
  """
  salmon index --transcripts $transcriptome --index salmon_index
  """
}
```

This process is now directed to copy all output files into a `results/`
directory. This saves having to specify the output directory in the script
definition each process, or a tedious `mv salmon_index/ results/` step. 
Nextflow also handles whether the directory already exists or if it
should be created. 

More information and other modes can be found on
[publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir).

> in the bash script you had to manually `mkdir -p "results`

## 2.1.3 Using Docker containers  

Nextflow recommends using containers to ensure reproducibility and portability
of your workflow. Containers package all the software and dependencies needed
for each tool into a self-contained environment. This means you don’t have to
manually install anything on your system, and your workflow will work
consistently across different systems — whether you're running it on your
local machine, a cluster, or in the cloud. Containers make it easier to share
your workflow with others and ensure it runs the same way every time, no matter
where it's executed.  

Nextflow supports
[multiple container runtimes](https://www.nextflow.io/docs/latest/container.html#).
In this workshop, we'll be demonstrating the value containers can bring to your
workflow by using Docker.

!!! info "Remember: different tools for different purposes"  

    In this workshop, we're using Docker to run containers. However, for some
    systems like HPC where you won't have administrative access to your environment,
    other options like Singularity/Apptainer will be more suitable.
    
    You don't have to write your own containers to run in your workflow. There are
    many container repositories out there. We highly recommend using 
    [Biocontainers](https://biocontainers.pro/registry) wherever possible.
    Biocontainers are pre-built and tested containers specifically for
    bioinformatics tools. They have a huge library and great community support. 
    
    You can find Biocontainers at the following repositories:  
    
    * [Biocontiners registry](https://biocontainers.pro/registry)
    * [Quay.io](https://quay.io/organization/biocontainers)
    * [DockerHub](https://hub.docker.com/r/biocontainers/biocontainers)
    * [Seqera containers](https://seqera.io/containers/)

Add the following container directive to the `INDEX` process, above
`publishDir`:  

In Nextflow, we can run a process using the [container](https://www.nextflow.io/docs/latest/process.html#container) directive, `container`.  
Add the following container directive to the `INDEX` process, above `publishDir`:  
```groovy title="main.nf"
process INDEX {

    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
    publishDir "results" mode: 'copy'

    input:
```

Usually, containers need to be downloaded using a command such as
`docker pull [image]`. All containers have been previously downloaded for the
workshop beforehand.  

??? tip "Use one container per process"
    
    Using single containers for each process in your workflow is considered best practices for the following reasons:

    - **Flexibility**: different processes require different tools (or versions). By using separate containers, you can easily tailor the container to the needs of each process without conflicts.
    - **Build and run efficiency**: Smaller, process-specific containers are faster to load and run compared to one large container that has unnecessary tools or dependencies for every process.
    - **Easier Maintenance**: it’s easier to update or modify one container for a specific process than to manage a large, complex container with many tools.
    - **Reproducibility**: reduces the risk of issues caused by software conflicts.

Before we can run the workflow, we need to tell Nextflow to run containers
using Docker. Nextflow requires [Docker](https://www.nextflow.io/docs/latest/container.html#docker)
to be installed on your system in order for this to work. Docker has been pre-installed on your Virtual Machine.  

We can tell Nextflow configuration to run containers with Docker by using a `nextflow.config` file.

Create a `nextflow.config` file in the same directory as `main.nf` and add the
following:

```groovy linenums="1" title="nextflow.config"
docker.enabled = true
```

You now have a complete process! 

## 2.1.4 Adding `params` and the workflow scope  

Now that you have written your first Nextflow process, we need to preapre it
for execution.  

You can think of Nextflow processes as similar to a function definition in R
or Python. We have defined what the process should do, but to actually run it,
we need to call the process within the workflow and pass in the inputs.

To run the process, we need to call it inside the `workflow{}` block, where
we control how data flows through the pipeline. To provide the input data we
need to define parameters. 

In the `00_index.sh` script, the file `data/ggal/transcriptome.fa` was passed
as the input into `salmon index`.  

We will pass in this file path as a `params`. Add the following to the top of your `main.nf` script:  

```groovy title="main.nf"

// pipeline input parameters
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
```

Recall that [parameters](https://www.nextflow.io/docs/latest/module.html#module-parameters)
are inputs and options that can be customised before the workflow is executed.

They allow you to control things like file paths, options for tools without
changing the process code itself.  

In our case, we could use `params.transcriptome` to provide a different
transcriptome file.

We define it in the `main.nf` script instead of in a command with `--` double
hyphen as the file (path) will not change.

[`$projectDir`](https://www.nextflow.io/docs/latest/script.html#configuration-implicit-variables)
is a configuration implicit variable that indicates the directory of the
`main.nf` script. 

Next, add the workflow scope at the bottom of you `main.nf` after the process:  

```groovy title="main.nf"
workflow {
    INDEX(params.transcriptome_file)
}
```

This will tell Nextflow to run the `INDEX` process with
`params.transcriptome_file` as input.

We are now ready to run our workflow!  

## 2.1.5 Running the workflow  

In the terminal, run the command:  

```bash
nextflow run main.nf
``` 

Your output should look something like:  

```console title="Output"
N E X T F L O W   ~  version 24.04.4

Launching `main.nf` [chaotic_jones] DSL2 - revision: 6597720332

executor >  local (1)
[de/fef8c4] INDEX | 1 of 1 ✔
```

Recall that the specifics of the output are randomly generated (i.e.
`[chaotic_jones]` and `[de/fef8c4]` in this example).

In this example, the output files for the `INDEX` process is output in
`work/26/c410b1...`.

> Inspect results folder?  

> results/salmon_index has a bunch of different files

!!! question "Exercise"

    Inspect the `.command.sh` file and compare it with `00_index.sh`. Note the
    similarities and differences.  

You have successfully run your first workflow!  

!!! abstract "Summary"

    In this step you have learned:  

    1. How to implement a simple process with input data  
    2. How to define parameters in your workflow scripts and the command
    line
    3. How to use a Docker container for a process  
    4. How to output files in a dedicated `publishDir`  

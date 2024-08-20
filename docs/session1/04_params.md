# Parameters

Parameters are constructs that can hold command line arguments.

They are useful because there will be many parameters, such as filenames and processing options, that you may want to decide or change at the time you execute your run. Simply, parameters allow you to do this without editing the script itself.

Here we're going to update the script with parameters to make it more flexible.

## Adding parameters

Parameters can be created by prefixing a parameter name with the parameters scope (e.g., `params.outdir`). They are accessible by processes and workflows and can be modified when you run your pipeline by adding a double hyphen (`--`) to the start of the parameter.

It is good practice to make a pipelines publishing directory flexible so outputs from a pipeline are not all put into the same place. As such, making the publishing directory a parameter would allow the user to rename it at the time of execution: 

```
publishDir params.outdir
```

The parameter `--outdir` could now be included in the run command (e.g., `nextflow run hello-world.nf --outdir results`).

!!!question "Exercise"

    Replace `results` with an `outdir` parameter in process block:

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="3"
        // Use echo to print 'Hello World!' and redirect to output.txt
        process SAYHELLO {
            publishDir params.outdir

            output: 
            path 'output.txt'
            
            script:
            """
            echo 'Hello World!' > output.txt
            """
        }
        ```

## Add a default value

No default value has been supplied to the pipeline.

Executing `hello-world.nf` without supplying a vale for `outdir` will throw an error.

```console
ERROR ~ Unexpected error while finalizing task 'SAYHELLO' - cause: Target path for directive publishDir cannot be null
```

Without a default value, a value must be supplied to `outdir` each execution.

```bash
nextflow run hello-world.nf --outdir new_results
```

!!!question "Exercise"

    Execute `hello-world.nf` with the `--outdir` parameters flag and a file name of your choice. View the output folder in your working directory.

    ???"Solution"

        ```bash
        nextflow run hello-world.nf --outdir new_results
        ```

        ```bash
        ls new_results
        ```

!!! abstract "Summary"

    In this step you have learned:  

    1. How to how to add a parameter to a pipeline 
    2. How to modify a parameter using the command line

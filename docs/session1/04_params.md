# Parameters

Parameters (`params` of short) are constructs that can hold command line arguments.

They are useful because there will be many parameters, such as filenames and processing options, that you may want to decide or change at the time you execute your run. Parameters allow you to do this without editing the script itself.

Parameters can be created by prefixing a parameter names with the params scope (e.g., `params.outdir`) and are accessible by processes and workflows.

Parameters can be modified when you run your pipeline by adding a double hyphen (`--`) to the start of the parameter name and including it in the run command (e.g., `nextflow run hello-world.nf --outdir results`).

## Adding parameters

Parameters are accessible in process blocks.

!!!question "Exercise"

    Replace `results` with an `outdir` parameter in process block:

    ??? "Solution"

        ```groovy title="hello-world.nf"
        process SAYHELLO {
            publishDir params.outdir

            output: 
            <truncated>
        ```

## Show changes

No default value has been supplied to the pipeline. Executing `hello-world.nf` without supplying a vale for `outdir` will throw an error.

```console
ERROR ~ Unexpected error while finalizing task 'SAYHELLO' - cause: Target path for directive publishDir cannot be null
```

Without a default value, a value must be supplied to `outdir` each execution.

```bash
nextflow run hello-world.nf --outdir param_results
```

!!!question "Exercise"

    Execute `hello-world.nf` with an outdir named `param_results`:

    ??? "Solution"

    ```bash
    nextflow run hello-world.nf --outdir param_results
    ```

    You should now see a new folder named `param_results` in your working directory.

!!! abstract "Summary"

    In this step you have learned:  

    1. How to how to add a parameter to a pipeline 
    2. How to modify a parameter using the command line

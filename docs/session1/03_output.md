# Outputs

Instead of printing "Hello World!" to the standard output it can be saved to a file. In a "real-world" pipeline, this is like having a command that specifies an output file as part of its normal syntax.

Here we're going to update the script and the output definition blocks to save the "Hello World!" output.

## Redirecting your output

The script block can be updated to redirect the `Hello World!` output to a file.

The `>` operator can be used for output redirection.

!!!question "Exercise"

    Redirect "Hello World!" to a file named `output.txt` in your script block

    ??? "Solution"

    ```groovy
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

## Outputs definitions

Outputs definitions in the process blocks typically require a qualifier and a variable name:

```groovy
output:
<output qualifier> <output name>
```

The **output qualifier** defines the type of data to be received. This information is used by Nextflow to apply the semantic rules associated with each qualifier, and handle it properly.
    
Common output qualifiers include `val` and `path`.

- `val`: Emit the variable with the specified name, e.g., `val <string>` 
- `path`: Emit a file produced by the process with the specified name,e.g., `path <file>`

If you set these wrong your pipeline will likely throw errors.

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#outputs) for a full list of output qualifiers.

The **output name** is a name given to the output variable. If a specific file is being produced it can be named in single quotes. 

```groovy title="hello-world.nf"
output:
path 'output.txt'
```

!!!question "Exercise"

    Execute the `hello-world.nf` pipeline again. Find the `output.txt` file in your work directory and verify that it contains the expected greeting.

!!! warning

    This example is brittle because we hardcoded the output filename in two separate places (the script and the output blocks). If we change one but not the other, the script will break.

## Publishing directory

Without a publishing strategy any files that are created by a process will only exist in the `work` directory. Realistically, you may want to capture a set of outputs and save them in a specific directory.

The `publishDir` directive can be used to specify where and how output files should be saved. For example:

```groovy
publishDir 'results'
```

By adding the following to a process, all output files would be saved in a new folder called `results` in the current working directory. 

!!!question "Exercise"

    Replace `debug true` with `publishDir 'outputs'` in your `SAYHELLO` process block. Execute your pipeline again. View the files in your working directory.

    ??? "Solution"

        ```groovy
        SAYHELLO {
            publishDir 'results'

            output:
            <truncated>
        ```

!!! abstract "Summary"

    In this step you have learned:  

    1. How to redirect outputs 
    2. How to use output definitions

# Use a parameter for naming the output file

Here we introduce `params` (short for 'parameters') as the construct that can hold command line arguments.

This is useful because there will be many parameters, such as filenames and processing options, that you may want to decide or change at the time you run your pipeline. Parameters allow you to do this without editing the script itself every time.

Parameters can be created by prefixing a parameter names with the params scope (e.g., `params.output_file`). 

Parameters can be modified when you run your pipeline by adding a double hyphen (`--`) to the start of the parameter name and including it in the run command (e.g., `nextflow run hello-world.nf --output_file results`).

## Adding an output file parameter

The `params` scope allows you to define parameters that will be accessible in the pipeline script.

Simply, parameters can be added to processes as inputs and outputs.

```groovy title="hello-world.nf"
output:
    path params.output_file
```

## Change the script to use the parameter

When including a parameter in a script block, a `$` must be used to treat it like a variable.

!!!question "Exercise"

    Add the `output_file` parameter to your `script` block:

    ??? "Solution"

        ```groovy title="hello-world.nf"
        script:
        """
        echo 'Hello World!' > $params.output_file
        """
        ```

## Show changes

```bash
nextflow run hello-world.nf --output_file 'output.txt'
```

The log output should start looking very familiar:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [evil_bose] DSL2 - revision: 6907ac9da2
executor >  local (1)
[46/e4ff05] process > SAYHELLO [100%] 1 of 1 ✔
```

Follow the same procedure as before to find the `output.txt` output file. If you want to convince yourself that the parameter is working as intended, feel free to repeat this step with a different output filename.

!!! warning

    If you forget to add the output filename parameter, you get a warning and the output file is called `null`. If you add it but don't give it a value, the output file is called `true`.

!!! abstract "Summary"

    In this step you have learned:  

    1. How to  
    2. How to 
    3. How to 

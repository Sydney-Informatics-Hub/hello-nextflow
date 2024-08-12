## 4. Use a command line parameter for naming the output file

Here we introduce `params` (short for 'parameters') as the construct that can hold command line arguments.

This is useful because there will be many parameters, such as filenames and processing options, that you may want to decide at the time you run the workflow. Parameters allow you to do this without editing the script itself every time.

Parameters can be created by prefixing a parameter names with the params scope (e.g., `params.output_file`). When including these in a script block, a `$` must be used to treat it like a variable.

Parameters can be modified when you run your workflow by adding a double hyphen (`--`) to the start of the parameter name and including it in the run command (e.g., `nextflow run hello-world --output_file results`).

## 1. Change the output declaration in the process to use a parameter

_Before:_

```groovy title="hello-world.nf" linenums="6"
output:
    path 'output.txt'
```

_After:_

```groovy title="hello-world.nf" linenums="6"
output:
    path params.output_file
```

#### 2. Change the process command to use the parameter too

_Before:_

```groovy title="hello-world.nf" linenums="9"
"""
echo 'Hello World!' > output.txt
"""
```

_After:_

```groovy title="hello-world.nf" linenums="9"
"""
echo 'Hello World!' > $params.output_file
"""
```

#### 3. Run the workflow again with the `--output_file` parameter

```bash
nextflow run hello-world.nf --output_file 'output.txt'
```

The log output should start looking very familiar:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [evil_bose] DSL2 - revision: 6907ac9da2
executor >  local (1)
[46/e4ff05] process > sayHello [100%] 1 of 1 ✔
```

Follow the same procedure as before to find the `output.txt` output file. If you want to convince yourself that the parameter is working as intended, feel free to repeat this step with a different output filename.

!!! warning

    If you forget to add the output filename parameter, you get a warning and the output file is called `null`. If you add it but don't give it a value, the output file is called `true`.

!!! tip

    Command-line arguments take one dash (-) for Nextflow options, two dashes (--) for workflow parameters.

### Takeaway

You know how to use a command line parameter to set the output filename.

### What's next?

Learn how to set a default value in case we leave out the parameter.

---

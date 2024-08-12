# Send the output to a file

Instead of printing "Hello World!" to the standard output it can be saved to a file.

In a real-world workflow, this is like having a command that specifies an output file as part of its normal syntax.

Both the script and the output definition blocks need to be updated.

!!! note

    Inputs and outputs in the process blocks typically require a qualifier and a variable name:

    ```
    <input/output qualifier> <input/output name>
    ```

    A definition consists of a qualifier and a name. The qualifier defines the type of data to be received. This information is used by Nextflow to apply the semantic rules associated with each qualifier, and handle it properly.
    
    Common qualifiers include `val` and `path`.

!!!question "Exercise"

    Write your output to a file named `output.txt`

    _Before:_

    ```groovy title="hello-world.nf" linenums="9"
    """
    echo 'Hello World!'
    """
    ```

    _After:_

    ```groovy title="hello-world.nf" linenums="9"
    """
    echo 'Hello World!' > output.txt
    """
    ```

## Change the output declaration in the process

_Before:_

```groovy title="hello-world.nf" linenums="6"
output:
    stdout
```

_After:_

```groovy title="hello-world.nf" linenums="6"
output:
    path 'output.txt'
```

## Run the workflow again

```bash
nextflow run hello-world.nf
```

The log output should be very similar to the first time your ran the workflow:

```console title="Output"
N E X T F L O W  ~  version 23.10.1
Launching `scripts/hello-world.nf` [disturbed_cajal] DSL2 - revision: 9512241567
executor >  local (1)
[ab/c61321] process > sayHello [100%] 1 of 1 ✔
```

Find the `work` directory in the file explorer. Find the `output.txt` output file and click on it to open it and verify that it contains the greeting as expected.

!!! warning

    This example is brittle because we hardcoded the output filename in two separate places (the script and the output blocks). If we change one but not the other, the script will break.

!!! abstract "Summary"

    In this step you have learned:  

    1. How to  
    2. How to 
    3. How to 

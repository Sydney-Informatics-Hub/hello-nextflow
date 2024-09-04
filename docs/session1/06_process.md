# Adding processes

Up until now you've been modifying a single step. However, pipelines generally consist of multiple steps where outputs from one step are used as inputs for the next.

Here you're going to step things up again and add another process to the pipeline.

## Translating text

The `tr` command is a UNIX command-line utility for **translating** or deleting characters. It supports a range of transformations including uppercase to lowercase, squeezing repeating characters, deleting specific characters, and basic find and replace. It can be used with UNIX pipes to support more complex translation. `tr` stands for translate. 

```bash
tr '[a-z]' '[A-Z]'`
```

## Piping commands

The pipe command in Linux, represented by the vertical bar symbol `|`, is an essential tool for command-line enthusiasts and professionals alike. The primary purpose of the pipe command is to connect the output of one command directly into the input of another:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]'
```

The contents of a file can be piped using the `cat` command:

```bash
cat output.txt | tr '[a-z]' '[A-Z]'
```

Like before, the output can be redirected to an output file:

```bash
cat output.txt | tr '[a-z]' '[A-Z]' > upper.txt
```

## `CONVERTTOUPPER`

The output of the `SAYHELLO` process is a text file called `output.txt`.

In the next step of the pipeline, you will add a new process named convert `CONVERTTOUPPER` that will convert all of the lower case letters in this file to a uppercase letters and save them as a new file.

The `CONVERTTOUPPER` process will follow the same structure as the `SAYHELLO` process:

```groovy
process CONVERTTOUPPER {
    publishDir params.outdir

    input:
    <input qualifier> <input name>

    output:
    <output qualifier> <output name>

    script:
    """
    <script>
    """
}
```

Using what you have learned in the previous sections you will now write a new process using the `tr` command from above.

!!!question "Exercise"

    Add new process named `CONVERTTOUPPER` that will take an input text file, convert all of the lowercase letters in the text file to uppercase letters, and save a new text file that contains the translated letters.

    ???Tip "Hint: `input:`"

        ```
        path input_file
        ```
        
        _Hint 1: The input is a file and requires the `path` qualifier._

        _Hint 2: The input name is `input_file`, however, you may call it something different._

    ???Tip "Hint: `output:`"

        The output

        ```
        path 'upper_output.txt'
        ```

        _Hint 1: The output is a file and requires the `path` qualifier._

        _Hint 2: The output name is hard coded as 'upper.txt', however, you may call it something different._

    ???Tip "Hint: `script:`"

        The script might look something like this:

        ```groovy
        cat $input_file | tr '[a-z]' '[A-Z]' > upper.txt
        ```

        _Hint 1: `input_file` must be the same as what was specified as the input name in the input block._

        _Hint 2: The output text file is named `upper.txt`

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="17-30"
        // Set default greeting
        params.greeting = 'Hello World!'

        // Set a default output directory
        params.outdir = 'results'

        // Use echo to print 'Hello World!' and redirect to output.txt
        process SAYHELLO {
            publishDir params.outdir

            input:
            val greeting

            output: 
            path 'output.txt'
            
            script:
            """
            echo '$greeting' > output.txt
            """
        }

        process CONVERTTOUPPER {
            publishDir params.outdir
            
            input:
                path input_file

            output:
                path 'upper.txt'

            script:
            """
            cat $input_file | tr '[a-z]' '[A-Z]' > upper.txt
            """
        }

        workflow {

            // Create a channel for inputs
            greeting_ch = Channel.of(params.greeting)

            // Emit a greeting
            SAYHELLO(greeting_ch)

        }
        ```

## Connecting processes

Outputs from one process can be used as inputs for another.

Outputs from a process can be accessed by adding `.out` to the end of a process name in the workflow definition:

```groovy
SAYHELLO.out
```

Outputs can then be used as an input for another process:

```groovy
CONVERTTOUPPER(SAYHELLO.out)
```

The same output could be used as inputs for multiple processes. 

!!!warning

    Adding `.out` to the end of a process name only works for single outputs. If there are multiple outputs the `emit` option must be used. See [additional options](https://www.nextflow.io/docs/latest/process.html#additional-options) for more information.

!!!question "Exercise"

    Add the `CONVERTTOUPPER` process to your workflow definition. Use the output from `SAYHELLO` as its input.

    ???solution

        ```groovy title="hello-world.nf" hl_lines="17-30"
        // Set default greeting
        params.greeting = 'Hello World!'

        // Set a default output directory
        params.outdir = 'results'

        // Use echo to print 'Hello World!' and redirect to output.txt
        process SAYHELLO {
            publishDir params.outdir

            input:
            val greeting

            output: 
            path 'output.txt'
            
            script:
            """
            echo '$greeting' > output.txt
            """
        }

        process CONVERTTOUPPER {
            publishDir params.outdir
            
            input:
                path input_file

            output:
                path 'upper.txt'

            script:
            """
            cat $input_file | tr '[a-z]' '[A-Z]' > upper.txt
            """
        }

        workflow {

            // Create a channel for inputs
            greeting_ch = Channel.of(params.greeting)

            // Emit a greeting
            SAYHELLO(greeting_ch)

            // Convert the greeting to uppercase
            CONVERTTOUPPER(SAYHELLO.out)

        }
        ```

Executing `hello-world.nf` will now show a second step:

```console
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [mighty_murdock] DSL2 - revision: 80e92a677c
executor >  local (2)
[ef/b99a2f] SAYHELLO (1)       [100%] 1 of 1 ✔
[cd/c8cf1b] CONVERTTOUPPER (1) [100%] 1 of 1 ✔
```

!!! abstract "Summary"

    In this step you have learned:  

    1. How to translate strings
    2. How add more processes to a script
    3. How to use outputs and inputs

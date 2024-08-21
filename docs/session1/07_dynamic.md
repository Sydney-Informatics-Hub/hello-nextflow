# Dynamic naming

Currently, the outputs of the `SAYHELLO` and `CONVERTTOUPPER` processes are being saved as `output.txt` and `upper.txt`, respectively.

In some situations this would be fine. However, to help identify the outputs we want our outputs to be dynamic.

Let's get tricky and name our output files dynamically.

### Dynamic outputs

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic string that references values defined in the input declaration block or in the script global context.

For example, the `SAYHELLO` input value `greeting` can be used to help name the output file.

```groovy hl_lines="5 8 12"
process SAYHELLO {
    publishDir 'results'

    input:
    val greeting

    output: 
    path "${greeting}.txt"
    
    script:
    """
    echo '$greeting' > ${greeting}.txt
    """
}
```

Curly brackets `{}` have been used to wrap `greeting` in the `output` and script definitions so it is interpreted as a variable as a part of a file name. 

There is an important difference between single-quoted (`'`) and double-quoted (`"`)strings. Double-quoted strings support variable interpolations while single-quoted strings do not.

!!!question "Exercise"

    Update the `SAYHELLO` and `CONVERTTOUPPER` process to use dynamic output names.

    !!!warning

        It's difficult to name a file with a space. Use a simple greeting, e.g., "Hello", when testing your pipeline.

    ???solution

        ```groovy title="hello-world.nf" hl_lines="9 13 25 29 36"
        // Use echo to print 'Hello World!' and redirect to output.txt
        process SAYHELLO {
            publishDir 'results'

            input:
            val greeting

            output: 
            path "${greeting}.txt"
            
            script:
            """
            echo '$greeting' > ${greeting}.txt
            """
        }

        // Use tr to convert lowercase letters to upper case letters and save as upper.txt
        process CONVERTTOUPPER {
            publishDir 'results'
            
            input:
                path input_file

            output:
                path "upper_${input_file}"

            script:
            """
            cat $input_file | tr '[a-z]' '[A-Z]' > upper_${input_file}
            """
        }

        workflow {

            // Set default greeting
            params.greeting = "Hello"

            // Create a channel for inputs
            greeting_ch = Channel.of(params.greeting)

            // Emit a greeting
            SAYHELLO(greeting_ch)

            // Convert the greeting to uppercase
            CONVERTTOUPPER(SAYHELLO.out)

        }
        ```

Let's execute our pipeline and view the changes to see if our outputs have been named dynamically.

```bash
nextflow run hello-world.nf --greeting Hello
```

While the output will look the same:

```console
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [mighty_murdock] DSL2 - revision: 80e92a677c
executor >  local (2)
[ef/b99a2f] SAYHELLO (1)       [100%] 1 of 1 ✔
[cd/c8cf1b] CONVERTTOUPPER (1) [100%] 1 of 1 ✔
```

We should now see some new files in our results folder:

```
Hello.txt       upper_Hello.txt
```

!!! abstract "Summary"

    In this step you have learned:  

    1. How to utilize dynamic naming
    2. How to use curly brackets (`{}`)
    3. How to use single (`'`) and double (`"`) quotes

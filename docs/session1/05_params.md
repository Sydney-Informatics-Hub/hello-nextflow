# Parameters

Parameters are constructs that can hold command line arguments.

They are useful because there will be many parameters, such as filenames and processing options, that we may want to change at the time we execute our run. Simply, parameters allow us to do this without editing the script itself.

Nextflow has multiple levels of configuration and as different levels may have conflicting settings, they are ranked in order of priority and some configuration can be overridden.

Parameters can be created by prefixing a parameter name with the parameters scope (e.g., `params.greeting`).

They are accessible by processes and workflows and can be modified when we run our pipeline by adding a double hyphen (`--`) to the start of the parameter.

Here we're going to update the script with parameters to make it more flexible.

## Greetings

Instead of hard coding "Hello World!" as a hard coded input, a parameter, with a default value can be created:

```groovy
params.greeting = "Hello World!"
```

The parameter can then be used in a channel factory, just like the hard coded string:

```groovy
greeting_ch = Channel.of(params.greeting)
```

The parameter can then be flexibly changed using a `--greeting flag` during execution:

```bash
nextflow run hello-world.nf --greeting 'Bonjour le monde!'
```

!!!question "Exercise"

    Update the workflow block with a default `greeting` parameter.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="3 4 7"
        workflow {

            // Set default greeting
            params.greeting = "Hello World!"

            // Create a channel for inputs
            greeting_ch = Channel.of(params.greeting)

            // Emit a greeting
            SAYHELLO(greeting_ch)
        }
        ```

!!!question "Exercise"

    Execute `hello-world.nf` again with the `--greeting` flag and a custom greeting.

    ???Solution

        ```bash
        nextflow run hello-world.nf --greeting 'Bonjour le monde!'
        ```

We can now check the output to see if our new `output.txt` file contains our new greeting.

!!! abstract "Summary"

    In this step we have learned:  

    1. How to how to add a parameter to a pipeline 
    2. How to modify a parameter using the command line

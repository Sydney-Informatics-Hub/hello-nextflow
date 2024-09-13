# Parameters

!!! info "Learning objectives"

    1. Implement pipeline parameters
    2. Understand the importance of parameters for flexible pipelines

Parameters are constructs that can hold command line arguments.

Here you're going to update the script with parameters to make it more flexible.

## Why are parameters useful?

Nextflow has multiple levels of configuration and, as different levels may have conflicting settings, they are ranked in order of priority and some configuration can be overridden.

Parameters are useful because they can be set with a default value in a script but can then be overwritten at runtime using a flag. Simply, parameters allow us to configure some aspect of a pipeline without editing the script itself.

Parameters can be created by prefixing a parameter name with the parameters scope (e.g., `params.greeting`) and are accessible by processes and workflows. They can be modified when you run your pipeline by adding a double hyphen (`--`) to the start of the parameter name (`--greeting`) and adding it to an execution command:

```bash
nextflow run hello-world.nf --greeting 'Hey'
```

## `--greeting`

Instead of hard coding 'Hello World!' as an input, a parameter, with a default value, can be created:

```groovy
params.greeting = 'Hello World!'
```

The parameter can then be used in a channel factory (just like the hard coded string):

```groovy
greeting_ch = Channel.of(params.greeting)
```

The parameter can then be flexibly changed using a `--greeting` flag in the run command:

```bash
nextflow run hello-world.nf --greeting 'Bonjour le monde!'
```

!!!question "Exercise"

    Update the `hello-world.nf` script to use a `greeting` parameter as an input. Define the default for the `greeting` parameter at the top of the script and give it the default value `'Hello World!'`.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="1 2 23"
        // Set default greeting
        params.greeting = 'Hello World!'

        // Use echo to print 'Hello World!' and redirect to output.txt
        process SAYHELLO {
            publishDir 'results'

            input:
            val greeting

            output:
            path 'output.txt'

            script:
            """
            echo '$greeting' > output.txt
            """
        }

        workflow {

            // Create a channel for inputs
            greeting_ch = Channel.of(params.greeting)

            // Emit a greeting
            SAYHELLO(greeting_ch)

        }
        ```

The `hello-world.nf` pipeline can now be executed with the `--greeting` flag and a custom greeting:

```bash
nextflow run hello-world.nf --greeting 'Bonjour le monde!'
```

## `--outdir`

It isn't very convenient to have the same output directory created every time you run your pipeline as the results are being overwritten.

Instead, a parameter can be used so you can change the publishing directory for every execution:

```groovy
publishDir params.outdir
```

A default value can be used for convenience as Nextflow will throw and error if `publishDir` is set to `null`.

```groovy
params.outdir = 'results'
```

However, you may consider having no default value here and letting the pipeline fail to prevent the accidental overwriting of results.

!!!question "Exercise"

    Update the `hello-world.nf` script to use an `outdir` parameter as the publishing directory. Define the default for the `outdir` parameter at the top of the script and give it the default value `'results'`.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="4 5 9"
        // Set default greeting
        params.greeting = 'Hello World!'

        // Set default output directory
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

        workflow {

            // Create a channel for inputs
            greeting_ch = Channel.of(params.greeting)

            // Emit a greeting
            SAYHELLO(greeting_ch)

        }
        ```

!!! abstract "Summary"

    In this step you have learned:

    1. How to how to add a parameter to a pipeline
    2. How to modify a parameter using the command line

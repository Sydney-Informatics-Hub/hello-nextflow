# Inputs

!!! info "Learning objectives"

    1. Describe Nextflow channel types
    2. Utlizie Nextflow process input blocks

So far, you've been emitting a greeting ('Hello World!') that has been hardcoded into the script block. In a more realistic situation, you might want to pass a variable input to your script, much like you pass files to command line tools for analysis.

Here you're going to to add some flexibility by introducing **channels** to your workflow and an **input definition** to your `SAYHELLO` process.

## Channels

In Nextflow, processes primarily communicate through **channels**.

Channels are created using channel factories.

There are numerous types of channel factories which can be utilized for creating different channel types and data types.

Importantly, there are two kinds of channels (**queue** channels and **value** channels) which behave differently.

**Queue channel**

- A non-blocking unidirectional first-in first-out queue connecting a producer process (i.e. outputting a value) to a consumer process, or an operators.
- Can be consumed only once.

**Value channel**

- Can be bound (i.e. assigned) with one and only one value.
- Can be consumed any number of times.

You're going to start by creating a channel that will contain your greeting with the `Channel.of()` channel factory.

!!!note

    You can build different kinds of channels depending on the shape of the input data.

### Channel.of()

The `Channel.of` method allows us to create a channel that emits the arguments provided to it. For example:

```groovy
ch_greeting = channel.of('Hello World!')
```

Would create a channel (named `ch_greeting`) that contains the string 'Hello World!'

Channels need to be created within the `workflow` definition.

!!!question "Exercise"

    Create a channel named `greeting_ch` with the 'Hello World!' greeting.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="2 3 4"
        workflow {

            // Create a channel for inputs
            greeting_ch = Channel.of('Hello world!')

            // Emit a greeting
            SAYHELLO()
        }
        ```

## Input definition blocks

Before `greeting_ch` can be passed to the `SAYHELLO` process as an input, you must first add an **input block** in the process definition.

The inputs in the input block, much like the output block, must have a qualifier and a name:

```
<input qualifier> <input name>
```

Input names can be treated like a variable, and while the name is arbitrary, it should be recognizable.

No quote marks are needed for variable inputs.

```
val greeting
```

!!!question "Exercise"

    Add an `input` block to the `SAYHELLO` process  with an input. Update the comment at the same time.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="1 4 5 6"
        // Use echo to print a string and redirect to output.txt
        process SAYHELLO {
            publishDir 'results'

            input:
            val greeting

            output:
            path 'output.txt'

            script:
            """
            echo 'Hello World!' > output.txt
            """
        }
        ```

The `SAYHELLO` process is now expecting an input value.

The `greeting_ch` channel can now be supplied to `SAYHELLO()` process within the workflow block:

```groovy
SAYHELLO(greeting_ch)
```

Without this, Nextflow will throw an error.

!!!question "Exercise"

    Add the `greeting_ch` as an input for the `SAYHELLO` process.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="7"
        workflow {

            // Create a channel for inputs
            greeting_ch = Channel.of('Hello world!')

            // Emit a greeting
            SAYHELLO(greeting_ch)
        }
        ```

The final piece is to update the `script` block to use the `input` value.

For an input to be treated like a variable in the script block, a `$` must be prepended to the input name:

```groovy
echo '$greeting' > output.txt
```

The `'` around `$greeting` are required to treat the greeting as a single string.

!!!question "Exercise"

    Update `hello-world.nf` to use the greeting input.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="13"
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
            greeting_ch = Channel.of('Hello world!')

            // Emit a greeting
            SAYHELLO(greeting_ch)
        }
        ```

!!!note

    **The number of inputs in the input block and the workflow must match!** If you had multiple inputs they would be listed across multiple lines in the process input block and listed inside the brackets in the workflow block.

    ???tip "Example"

        ```groovy title="example.nf"
        process MYFUNCTION {
            debug true

            input:
            val input_1
            val input_2

            output:
            stdout

            script:
            """
            echo $input_1 $input_2
            """
        }

        workflow {
            MYFUNCTION('Hello', 'World!')
        }
        ```

**Yes! Your pipeline now uses an input channel!**

!!! abstract "Summary"

    In this step you have learned:

    1. How to use Channel factories
    2. How to how to add process inputs

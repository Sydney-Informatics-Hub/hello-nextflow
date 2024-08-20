# Inputs

So far, we've been emitting a greeting hardcoded into the process command.

Now we're going to add some flexibility by introducing channels.

## Channels

In Nextflow, processes primarily communicate through channels.

Channels are created used channel factories methods.

There are several types of channel factories which can be utilized for creating different channel types for different data types.

Importantly, there are two kinds of channels (**queue** channels and **value** channels) which behave differently.

**Queue channel**

- A non-blocking unidirectional first-in first-out queue connecting a producer process (i.e. outputting a value) to a consumer process, or an operators
Can be consumed only once

**Value channel**

- Can be bound (i.e. assigned) with one and only one value
- Can be consumed any number of times

We're going to start by creating a channel with the `Channel.of()` channel factory.

!!!note

    We can build different kinds of channels depending on the shape of the input data.

## Channel.of()

The `Channel.of` method allows us to create a channel that emits the arguments provided to it. For example:

```groovy
ch = channel.of( 1, 3, 5, 7 )
```

Would create a channel (named `ch`) that contains four values (1, 3, 5, and 7).

Channels need to be created within the `workflow` block.

## Input

The `input` definition in the process block, much like the `output` definition, must have a qualifier and a name.

```
<input qualifier> <input name>
```

Input names can be treated like a variable, and while the name is arbitrary, it should be recognizable.

!!!question "Exercise"

    Add a path `input` definition to the `SAYHELLO` process block. Update the comment at the same time.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="1 5 6"
        // Use echo to print a string and redirect to output.txt
        process SAYHELLO {
            publishDir params.outdir

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

The `SAYHELLO` process is now expecting an input.

Without this, Nextflow will throw an error.

The `Channel.of` channel factory can be used to create a channel containing a greeting.

```groovy
greeting_ch = Channel.of('Hello world!')
```

The `greeting_ch` channel can then be supplied to `SAYHELLO()` process within the workflow block.


```groovy
SAYHELLO(greeting_ch)
```

Note how the channel is being used as the input for `SAYHELLO`.

!!!question "Exercise"

    Create a channel named `greeting_ch` with the "Hello World!" greeting and use it as an input for the `SAYHELLO` process.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="3 4 7"
        workflow {

            // Create a channel for inputs
            greeting_ch = Channel.of('Hello world!')

            // Emit a greeting
            SAYHELLO(greeting_ch)
        }
        ```

The final piece is to update the `scrip`t block to use the `input` value.

For an input to be treated like a variable in the script block, a `$` must be prepended to the input name.

```groovy
echo '$greeting' > output.txt
```

!!!question "Exercise"

    Update `hello-world.nf` to use the greeting input.

    ???Solution

        ```groovy title="hello-world.nf" hl_lines="13"
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
            greeting_ch = Channel.of('Hello world!')

            // Emit a greeting
            SAYHELLO(greeting_ch)
        }
        ```

!!!note

    The number of inputs must match! If we had multiple inputs they would be listed across multiple lines in the process input definition and listed inside the brackets in the workflow block.

    ???tip "example.nf"

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
            MYFUNCTION("Hello", "World!")
        }
        ```

## Parameter inputs

A parameter can also be used as an input.

Nextflow has multiple levels of configuration and as different levels may have conflicting settings, they are ranked in order of priority and some configuration can be overridden.

Importantly, parameters that are hard coded into scripts are overwritten by parameters supplied on the command line.

Instead of hard coding "Hello World!" as an input, a parameter, with a default value can be created:

```groovy
params.greeting = "Hello World!"
```

The parameter can then be used in a channel factory, just like the hard coded string:

```
greeting_ch = Channel.of(params.greeting)
```

!!!question "Exercise"

    Update the workflow block with a default `greeting` parameter and execute it with the greeting `Bonjour le monde!`

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

    Execute hello-world.nf again with a custom greeting.

    ???Solution

        ```bash
        nextflow run hello-world.nf --greeting 'Bonjour le monde!'
        ```

!!! abstract "Summary"

    In this step we have learned:  

    1. How to how to add process inputs
    2. How to use Channel factories 
    3. How to parameter inputs

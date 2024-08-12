# Interpret `hello-world.nf`

Nextflow scripts are built up of multiple parts.

A **process** is the basic processing primitive to execute a user script.

The process definition starts with the keyword `process`, followed by process name and finally the process body delimited by curly braces. The process body must contain a script block which represents the command or, more generally, a script that is executed by it.

A process may contain any of the following definition blocks: directives, inputs, outputs, when clauses, and of course, the script.

A **workflow** is a composition of processes and dataflow logic.

The workflow definition starts with the keyword `workflow`, followed by an optional name, and finally the workflow body delimited by curly braces.

Processes are connected through asynchronous first-in, first-out (FIFO) queues, called **channels**. The interaction between processes, and ultimately the workflow execution flow itself, are defined by the process input and output declarations.

Let's open the `hello-world.nf` script and look at how it's structured.

## Commenting your pipeline

The first block of code describes a **process** called `SAYHELLO` that writes its **output** to `stdout`:

```groovy title="hello-world.nf"
process SAYHELLO {
    debug true

    output:
        stdout

    """
    echo 'Hello World!'
    """
}
```

The `debug` directive is used so the `stdout` will print to the terminal.

The second block of code describes the **workflow** itself, which consists of one call to the `SAYHELLO` process.

```groovy title="hello-world.nf"
workflow {
    SAYHELLO()
}
```

## Comment your pipeline

It is always worthwhile to comment your pipelines so you and others can understand the code.

A single line comment can be added by prepending it with two forward slash.

`// This is my comment` 

!!!question "Exercise"

    Add comment to document what each your process and workflow blocks are doing:

    ```groovy title="hello-world.nf" linenums="1"
    /*
    * Use echo to print 'Hello World!' to standard out
    */
    process sayHello {
    ```

#### 3. Add an in-line comment above the process call

```groovy title="hello-world.nf" linenums="14"
workflow {

    // emit a greeting
    sayHello()
}
```

!!! abstract "Summary"

    In this step you have learned:  

    1. How to  
    2. How to 
    3. How to 

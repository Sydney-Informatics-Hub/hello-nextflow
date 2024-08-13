# Interpret `hello-world.nf`

Nextflow scripts are built up of multiple parts.

A **process** is the basic processing primitive to execute a user script.

The process definition starts with the keyword `process`, followed by process name, and finally the process body delimited by curly braces. The process body must contain a script block which represents the command or, more generally, a script that is executed by it.

A process may contain any of the following definition blocks: `directives`, `inputs`, `outputs`, `when` clauses, and of course, the `script`.

A **workflow** is a composition of processes and dataflow logic.

The workflow definition starts with the keyword `workflow`, followed by an optional name, and finally the workflow body delimited by curly braces.

Processes are connected through queues, called **channels**. The interaction between processes, and ultimately the workflow execution flow itself, are defined by the process input and output declarations in each process.

Let's open the `hello-world.nf` script again and review its structure.

## Commenting your pipeline

The first block of code describes a **process** called `SAYHELLO` with three definitions.

- **debug true**: a directive that will print the output to your console
- **output**: directing `script` outputs to be printed to `stdout`
- **script**: the `echo 'Hello World!'` command

Using `debug true` and `stdout` in combination will cause "Hello World!" to be printed to your terminal.

```groovy title="hello-world.nf"
process SAYHELLO {
    debug true

    output:
        stdout

    script:
    """
    echo 'Hello World!'
    """
}
```

The second block of code describes the **workflow** itself, which consists of one call to the `SAYHELLO` process.

```groovy title="hello-world.nf"
workflow {
    SAYHELLO()
}
```

It is always worthwhile to comment your pipelines so you and others can understand the code.

A single line comment can be added by prepending it with two forward slash (`//`).

```groovy
// This is my comment` 
```

Multi-line comments can be added using the following format:

```groovy
/*
 * Use echo to print 'Hello World!' to standard out
 */
```

It is up to you, as a developer, to choose how and where to comment your code.

!!!question "Exercise"

    Add a comment to your pipeline to describe what the **process** block is doing:

    ??? "Solution"

        Your solution may look something like this:

        ```groovy title="hello-world.nf"
        /*
         * Use echo to print 'Hello World!' to standard out
         */
        process SAYHELLO {
        ```

!!! abstract "Summary"

    In this step you have learned:  

    2. How to interpret `hello-world.nf`
    3. How to add comments to your pipelines

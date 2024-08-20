# Your first pipeline

Nextflow is a workflow orchestration engine and domain-specific language (DSL) that makes it easy to write data-intensive computational workflows.

## `hello-world.nf`

In Nextflow, **process** is the basic processing primitive to execute a user script.

The process definition starts with the keyword `process`, followed by process name, and finally the process body delimited by curly braces. The process body must contain a script block which represents the command or, more generally, a script that is executed by it.

A process may contain any of the following definition blocks: `directives`, `inputs`, `outputs`, `when` clauses, and of course, the `script`.

A **workflow** is a composition of processes and dataflow logic.

The workflow definition starts with the keyword `workflow`, followed by an optional name, and finally the workflow body delimited by curly braces.

Processes are connected through queues, called **channels**. The interaction between processes, and ultimately the workflow execution flow itself, are defined by the process input and output declarations in each process.

Let's review the structure of  `hello-world.nf`.

```groovy title="hello-world.nf" linenums="1"
process SAYHELLO {
    debug true

    output: 
    stdout
    
    script:
    """
    echo 'Hello World!'
    """
}

workflow {
    SAYHELLO()
}
```

The first block of code (lines 1-11) describes a **process** called `SAYHELLO` with three definitions.

- **debug**: a directive that, when true, will print the output to your console
- **output**: directing `script` outputs to be printed to `stdout` (standard output)
- **script**: the `echo 'Hello World!'` command

The second block of code (13-15) lines describes the **workflow** itself, which consists of one call to the `SAYHELLO` process.

!!!note

    Using `debug true` and `stdout` in combination will cause "Hello World!" to be printed to your terminal.

## Commenting your code

It is worthwhile to comment your code so you, and others, can easily understand your code.

In Nextflow, a single line comment can be added by prepending it with two forward slash (`//`):

```groovy
// This is my comment` 
```

Similarly, multi-line comments can be added using the following format:

```groovy
/*
 * Use echo to print 'Hello World!' to standard out
 */
```

As a developer you can to choose how and where to comment your code.

!!!question "Exercise"

    Add a comment to your pipeline to describe what the **process** block is doing:

    ??? "Solution"

        Your solution may look something like this:

        ```groovy title="hello-world.nf"
        /*
         * Use echo to print 'Hello World!' to standard out
         */
        process SAYHELLO {
        <truncated>
        ```

        Or this:

        ```groovy title="hello-world.nf"
        // Use echo to print 'Hello World!' to standard out
        process SAYHELLO {
        <truncated>
        ```

---

## Run `hello-world.nf`

The `nextflow run` command is used to execute pipelines.

```bash
nextflow run <pipeline.nf>
```

When a pipeline is stored locally you need to supply the full path to the script. However, if the pipeline has been submitted to GitHub (and you have an internet connection) you can execute it without a local copy. For example, the `hello` repository hosted on the `nextflow-io` GitHub account:

```bash
nextflow run nextflow-io/hello
```

!!!question "Exercise"

    Use ` nextflow run` command to execute `hello-world.nf`:

        ```bash
        nextflow run hello-world.nf
        ```

**Congratulations! You have just ran your first pipeline!**

Your console should look something like this:

```console linenums="1"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [mighty_murdock] DSL2 - revision: 80e92a677c
executor >  local (1) // (3)!
[4e/6ba912] process > SAYHELLO [100%] 1 of 1 ✔
Hello World!
```

**Line:**

1. The version of Nextflow that was executed 
2. The script and version names
3. The executor used (in the above case: local)
4. The first process is executed once, which means there is one task. The line starts with a unique hexadecimal value, and ends with the task completion information
5. The result string from stdout is printed

When a Nextflow pipeline is executed, a `work` directory is created. Processes are executed in isolated task directories. Each task uses a unique directory based on its [hash](https://www.nextflow.io/docs/latest/cache-and-resume.html#task-hash) (e.g., `4e/6ba912`) within the work directory.

When a task is created, Nextflow stages the task input files, script, and other helper files into the task directory. The task writes any output files to this directory during its execution, and Nextflow uses these output files for downstream tasks and/or publishing.

!!!note

    You can execute `tree work` to view the work directory structure.

!!! warning

    Your work directory might not have the same hash as the one shown above.

A series of files log files and any outputs are created by each task in the work directory:

-   **`.command.begin`**: Metadata related to the beginning of the execution of the process task
-   **`.command.err`**: Error messages (stderr) emitted by the process task
-   **`.command.log`**: Complete log output emitted by the process task
-   **`.command.out`**: Regular output (stdout) by the process task
-   **`.command.sh`**: The command that was run by the process task call
-   **`.exitcode`**: The exit code resulting from the command

These files are created by Nextflow to manage the execution of your pipeline. While these file are not required now, you may need to interrogate them to troubleshoot issues later.

As these are dot files you may need to use `ls -la` to view them.

!!!question "Exercise"

    Browse the `work` directory and view the `.command.sh` file.

    ??? "Solution"

        _Note: Your hash will be different to the example shown below_

        ```bash
        cat work/4e/6ba9138vhsbcbsc83bcka/.command.sh
        ```

!!! abstract "Summary"

    In this step you have learned:  

    1. How to create a Nextflow pipeline
    2. How to interpret `hello-world.nf`
    3. How to add comments to your pipelines 
    4. How to `run` a Nextflow pipeline
    5. How to view log files create by Nextflow

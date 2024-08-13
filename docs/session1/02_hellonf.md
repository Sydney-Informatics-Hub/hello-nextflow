# Your first pipeline

Nextflow is a workflow orchestration engine and domain-specific language (DSL) that makes it easy to write data-intensive computational workflows.

It is designed around the idea that the Linux platform is the lingua franca of data science. Linux provides many simple but powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations.

Here we're going to create and run a script (named `hello-world.nf`) that will print 'Hello World!'.

## `hello-world.nf`

Nextflow provides a robust command line interface for the management and execution pipelines. You can view a full list of nextflow options and commands using the `-h` option:

```bash
nextflow -h
```

Nextflow scripts are saved with the `.nf` suffix. For example, `myscript.nf`.

The `nextflow run` command is used to execute a local or remote pipeline.

```bash
nextflow run <pipeline.nf>
```

!!!question "Exercise"

    Copy the following text and add it to a file named `hello-world.nf`, save it, and execute it using the ` nextflow run` command:

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

    workflow {
        SAYHELLO()
    }
    ```

    ??? "Solution"

        ```bash
        nextflow run hello-world.nf
        ```

## What just happened?

Congratulations, you ran your first Nextflow pipeline!

You console should look something like this:

```console title="Output" linenums="1"
N E X T F L O W  ~  version 23.10.1 // (1)!
Launching `hello-world.nf` [mighty_murdock] DSL2 - revision: 80e92a677c // (2)!
executor >  local (1) // (3)!
[4e/6ba912] process > SAYHELLO [100%] 1 of 1 ✔ // (4)!
Hello World! // (5)!
```

1. The version of Nextflow that was executed
2. The script and version names
3. The executor used (in the above case: local)
4. The first process is executed once, which means there is one task. The line starts with a unique hexadecimal value, and ends with the task completion information
5. The result string from stdout is printed

When a Nextflow pipeline is run a `work` directory that stores various files is created.

Each task uses a unique directory based on its [hash](https://www.nextflow.io/docs/latest/cache-and-resume.html#task-hash) (e.g., `4e/6ba912`) within the work directory.

When a task is created, Nextflow stages the task input files, script, and other helper files into the task directory. The task writes any output files to this directory during its execution, and Nextflow uses these output files for downstream tasks and/or publishing.

!!! warning

    Your work directory won't necessarily have the same hash as the one shown above.

You can browse the `work` directory to find the log files and any outputs created by the task.

You should find the following files:

-   **`.command.begin`**: Metadata related to the beginning of the execution of the process task
-   **`.command.err`**: Error messages (stderr) emitted by the process task
-   **`.command.log`**: Complete log output emitted by the process task
-   **`.command.out`**: Regular output (stdout) by the process task
-   **`.command.sh`**: The command that was run by the process task call
-   **`.exitcode`**: The exit code resulting from the command

While these file are not required now, you may need to interrogate them to troubleshoot issues later.

!!!question "Exercise"

    Browse the work directory and view the `.command.sh` file.

    ??? "Solution"

        _Note: Your hash will be different to the example shown below_

        ```bash
        cat work/4e/6ba9138vhsbcbsc83bcka/.command.sh
        ```

!!! tip

    Some of the specifics will be different in your log output. For example, here `[mighty_murdock]` and `[4e/6ba912]` are randomly generated names, so those will be different every time.

!!! abstract "Summary"

    In this step you have learned:  

    1. How to create a Nextflow pipeline 
    2. How to run a Nextflow pipeline
    3. How to view log files create by Nextflow

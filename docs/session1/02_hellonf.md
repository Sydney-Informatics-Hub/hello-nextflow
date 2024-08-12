# `hello-world.nf`

`hello-world.nf`

## Your first pipeline

Copy the following text and add it to a file named `hello-world.nf` in your text editor:

```groovy title="hello-world.nf"
process SAYHELLO {
    debug true

    output: 
        stdout
    
    """
    echo 'Hello World!'
    """
}

workflow {
    SAYHELLO()
}
```

## Nextflow run

Execute your new pipeline with the `nextflow run` command:

```bash
nextflow run hello-world.nf
```

You console should look something like this:

```console title="Output" linenums="1"
N E X T F L O W  ~  version 23.10.1
Launching `hello-world.nf` [mighty_murdock] DSL2 - revision: 80e92a677c
executor >  local (1)
[4e/6ba912] process > sayHello [100%] 1 of 1 ✔
Hello World!
```

Congratulations, you ran your first Nextflow pipeline!

The most important output here is line 4, which reports that the `sayHello` process was successfully executed `1 of 1`.

When a Nextflow pipeline is run a `work` directory that stores various files is created.

Each task uses a unique directory based on its [hash](https://www.nextflow.io/docs/latest/cache-and-resume.html#task-hash) (e.g., `4e/6ba912`) within the work directory.

When a task is created, Nextflow stages the task input files, script, and other helper files into the task directory. The task writes any output files to this directory during its execution, and Nextflow uses these output files for downstream tasks and/or publishing.

!!! warning

    Your work directory won't necessarily have the same hash as the one shown above.

Browse the `work` directory in the file explorer to find the log files and any outputs created by the task. You should find the following files:

-   **`.command.begin`**: Metadata related to the beginning of the execution of the process task
-   **`.command.err`**: Error messages (stderr) emitted by the process task
-   **`.command.log`**: Complete log output emitted by the process task
-   **`.command.out`**: Regular output (stdout) by the process task
-   **`.command.sh`**: The command that was run by the process task call
-   **`.exitcode`**: The exit code resulting from the command

In this case, look for your output in the `.command.out` file.

!!! tip

    Some of the specifics will be different in your log output. For example, here `[mighty_murdock]` and `[4e/6ba912]` are randomly generated names, so those will be different every time.

!!! abstract "Summary"

    In this step you have learned:  

    1. How to  
    2. How to 
    3. How to 

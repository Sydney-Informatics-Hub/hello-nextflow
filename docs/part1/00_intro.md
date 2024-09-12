# Welcome

In **Part 1** of this workshop, you will create a toy multi-step Nextflow workflow.

We will start by exploring common bash commands. Next, we will turn these into a small single step Nextflow pipeline that will print a greeting to our terminal. In a series of exercises, we will then iterate on this pipeline to make it more flexible using outputs, inputs, and parameters. Finally, we will add a second step to the pipeline to turn our greeting into uppercase letters and name pipeline outputs dynamically.

During **Part 2**, the skills and concepts you have learned will be applied in a more realistic scenario.

## Create a work directory

It is good practice to organize projects into their own folders to make it easier to track and replicate experiments over time.

!!!question "Exercise"

    Create a new directory for all of today’s activities and move into it:

    ```bash
    mkdir ~/part1 && cd $_
    ```

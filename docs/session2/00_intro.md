# 2.0 Introduction  

> Explain the data  

To demonstrate a real-world biomedical scenario, you will implement a proof of concept RNA-Seq workflow which:  
1. Indexes a transcriptome file  
2. Performs quality controls  
3. Performs quantification  
4. Creates a MultiQC report  

This will be done by iteratively building the workflow by adding to a single `main.nf` script, and a `nextflow.config` file. 

The workflow will make use of commonly used bioinformatics tools:  
1. [Salmon](https://combine-lab.github.io/salmon/) is a tool for quantifying molecules known as transcripts through RNA-seq data.  
2. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a tool for quality analysis of high throughput sequence data. You can think of it as a way to assess the quality of your data.  
3. [MultiQC](https://multiqc.info/) searches a given directory for analysis logs and compiles an HTML report for easy viewing. It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.  

We will not explore how these tools work much further. These common tools will be used to teach you how to create a Nextflow workflow that chains together the inputs/outputs of files.  

## 2.0.1 Workshop structure  

In small blocks, you will iteratively build a workflow that runs:  
1. A single process with input parameters  
2. All the tools, start-to-finish, for a single sample  
3. Optimally for multiple samples/files  

Each added "block" will involve:  
1. Adding a new Nextflow process/bioinformatics tool to the workflow  
2. Adding channels to process input files with groovy  
3. Applying a (new) concept/feature/trick  

Concepts that have been previously covered will appear as prompts to be completed independently.  

## 2.0.2 Pipeline structure and design 

High level description (and diagram) of the workflow we're building, why its structured this way, inputs and outputs, etc. 

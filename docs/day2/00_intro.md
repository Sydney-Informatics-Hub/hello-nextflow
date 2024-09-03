# 2.0 Introduction  

In this workshop, we will create a multi-sample Nextflow workflow for
preparing RNAseq data. We will build the workflow, step-by-step, by
converting a series of provided bash scripts into small workflow components.  

Along the way, you will encounter Nextflow concepts (from yesterday, and some
new) and best practices for developing your own pipeline.  

## 2.0.1 From bash scripts to scalable workflows  

Imagine you are a bioinformatician in a busy research lab. Your team will be
receiving a large batch of samples that need to be processed through a series
of analysis steps.  

You have inherited a set of bash scripts from a former colleague, which were
used to process a handful of samples manually. These scripts are robust
and well-tested, but they were not designed with scalability in mind.  

As more samples come in, running these scripts one by one is becoming
increasingly tedious and error-prone.  

**You need a way to automate this process, ensuring consistency and efficiency
across many samples.**  

You decide to use Nextflow.  

## 2.0.2 RNAseq: data and tools  

RNAseq is used to study gene expression and has many applications across
biomedicine, agriculture and evolutionary studies.  

**Data**  

The data we will use includes:
- `*.fastq`: Paired-end RNAseq reads from three different samples (gut, liver,
lung)
- `transcriptome.fa`  

**Tools**  

We will be implementing and integrating three commonly used bioinformatics
tools:  
1. [Salmon](https://combine-lab.github.io/salmon/) is a tool for quantifying molecules known as transcripts through RNA-seq data.  
2. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a tool for quality analysis of high throughput sequence data. You can think of it as a way to assess the quality of your data.  
3. [MultiQC](https://multiqc.info/) searches a given directory for analysis logs and compiles an HTML report for easy viewing. It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.  

We will not explore how the data and tools work further, and focus on how
they should be implemented in a Nextflow workflow.  

## 2.0.3 Pipeline structure and design 

![](img/0.excalidraw.png)

High level description (and diagram) of the workflow we're building, why its structured this way, inputs and outputs, etc. 

Additionally, will need a clear explanation of how we plan to work through the materials. Explaining the .nf files, building workflow in one file, how we will be running software (containers). 

We will be following the `process` structure from Section 1.2:  

```groovy
process < name > {
  [ directives ]

  input:
    < process inputs >

  output:
    < process outputs >

  script:
  """
  < script to be executed >
  """
}
```

## 2.0.4 Advanced: `template-nf`  

Low barrier to developing workflows whilst adhering to best practices.

The workshop draws on _some_ of the concepts here.  

Should also provide a callout block that points to template-nf for more complex and modular workflow design.

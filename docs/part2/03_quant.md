# 2.3 Multiple inputs  

!!! note "Learning objectives"

    1. Implement a process with multiple inputs. 
    2. Practice implementing a Nextflow process from a bash script.  
    3. Practice chaining Nextflow processes with channels.  
    4. Understand the components of a process such as `input`, `output`,
    `script`, `directive`, and the `workflow` scope.  

> Update learning objective. Inputs outputs previously explained, adjust this one

In this section you will implement the `QUANTIFICATION` process. This section
contains more self-directed exercises for you to practice applying the
concepts introduced in the previous sections.

> Add more detail here, less focus on quantification, more on nextflow components they will be practicing

![](img/3.excalidraw.png)

You will implement the `02_quant.sh` bash script as a Nexflow process called
`QUANTIFICATION`.  

Open the bash script `02_quant.sh`.  

```bash title="02_quant.sh"
SAMPLE_ID=gut
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"

salmon quant \
    --libType=U \
    -i results/salmon_index \
    -1 ${READS_1} \
    -2 ${READS_2} \
    -o results/${SAMPLE_ID}
```

- `--libType=U` is a required argument. Leave this as is for the script definition.  
- `-i results/salmon_index` is the directory output by the `INDEX` process. **This is the output of the `INDEX` process.**
- `-1` and `-2` are flags for the respective paired reads (`.fq`). **These are output by the `reads_in` channel.**
- `-o` indicates that the command will output files into a directory called `results/gut`

## 2.3.1 Implementing the process  

### 1. Process directives  

Here is the empty `process` template with the `container` and `publishDir`
directives we'll be using to get you started. Add this after where you 
defined the `FASTQC` process.  

```groovy title="main.nf"
process QUANTIFICATION {
  
  container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
  publishDir "results", mode: 'copy'

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

### 2. Define the process `script`  

Update the `script` definition:

```groovy title="main.nf" hl_lines="14"
process QUANTIFICATION {

  container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
  publishDir "results", mode: 'copy'

  input:
    < process inputs >

  output:
    < process outputs >

  script:
  """
  salmon quant --libType=U -i $salmon_index -1 $reads_1 -2 $reads_2 -o $sample_id
  """
}
```

The `--libType=U` argument stays the same.

> Explain other variables: converting hardcoded paths and bash variables to
Nextflow `process` variables, dropping results/ due to publishDir results 

### 3. Define the process `output`

The `output` is a directory of `$sample_id`. In this case, it will be a
directory called `gut/`. Replace `< process outputs >` with the following:  

```groovy title="main.nf" hl_lines="10"
process QUANTIFICATION {

  container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
  publishDir "results", mode: 'copy'

  input:
    < process inputs >

  output:
  path "$sample_id"

  script:
  """
  salmon quant --libType=U -i $salmon_index -1 $reads_1 -2 $reads_2 -o $sample_id
  """
}
```

### 4. Define the process `input`  

In this step we will define the process inputs. Based on the bash script, we
have four inputs:  

- `$salmon_index`
- `$sample_id`
- `$reads_1`
- `$reads_2`

These should look familiar! 

The `$salmon_index` was output by the `INDEX` process, and `$sample_id`,
`$reads_1`, `$reads_2` are output by our `reads_pairs_ch`. We will see how to
chain these when we work on the `workflow` scope below.

First, add the input definition for `$salmon_index`. Recall that we use the
`path` qualifier as it is a directory:  

```groovy title="main.nf" hl_lines="7"  
process QUANTIFICATION {

  container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
  publishDir "results", mode: 'copy'

  input:
  path salmon_index

  output:
  path "$sample_id"

  script:
  """
  salmon quant --libType=U -i $salmon_index -1 $reads_1 -2 $reads_2 -o $sample_id
  """
}
```

Secondly, add the tuple input:  

> Needs more explanation i.e. use of val vs path

```groovy title="main.nf" hl_lines="8"
process QUANTIFICATION {

  container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
  publishDir "results", mode: 'copy'

  input:
  path salmon_index
  tuple val(sample_id), path(reads_1), path(reads_2)

  output:
  path "$sample_id"

  script:
  """
  salmon quant --libType=U -i $salmon_index -1 $reads_1 -2 $reads_2 -o $sample_id
  """
}
```

You have just defined a process with multiple inputs!  

### 5. Call the process in the `workflow` scope  

> More explanation and guidance needed. Guide through making the input channel from the indexing step.

Recall that the inputs for the `QUANTIFICATION` process are emitted by the
`reads_in` channel and the output of the `INDEX` process. The `reads_in`
is ready to be called by the `QUANTIFICATION` process, but we need to first
prepare an channel for the index files.  

Add the following 

```groovy title="main.nf" hl_lines="12-14"
// Define the workflow
workflow {

    // Run the index step with the transcriptome parameter
    INDEX(params.transcriptome_file)

    // Define the fastqc input channel
    reads_in = Channel.fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }

    // Define the quantification channel for the index files
    transcriptome_index_in = INDEX.out[0]

}
```

> connect back to part 1 https://training.nextflow.io/basic_training/processes/#output-definitions  

Call the `QUANTIFICATION` process in the workflow scope and add the inputs:  

```groovy title="main.nf" hl_lines="15-17"
// Define the workflow
workflow {

    // Run the index step with the transcriptome parameter
    INDEX(params.transcriptome_file)

    // Define the fastqc input channel
    reads_in = Channel.fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }

    // Define the quantification channel for the index files
    transcriptome_index_in = INDEX.out[0]

    // Run the quantification step with the index and reads_in channels
    QUANTIFICATION(transcriptome_index_in, reads_in)

}
```

We need to pass two arguments to `QUANTIFICATION` as there are two inputs in
the `process` definition. 

> Add note about tuples being a single input?

Run the workflow:  

```bash
nextflow run main.nf -resume
```

```console title="Output"
Launching `main.nf` [shrivelled_cuvier] DSL2 - revision: 4781bf6c41

executor >  local (1)
[de/fef8c4] INDEX              | 1 of 1, cached: 1 ✔
[bb/32a3aa] FASTQC (1)         | 1 of 1, cached: 1 ✔
[a9/000f36] QUANTIFICATION (1) | 1 of 1 ✔

```

You now have a `results/gut` folder, with an assortment of files and
directories.

!!! abstract "Summary"

    In this step you have learned:

        1. How to          
        1. How to 
        1. How to 

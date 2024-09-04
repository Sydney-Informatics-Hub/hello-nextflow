# 2.3 Multiple inputs  

!!! note "Learning objectives"

    1. Implement a process with multiple inputs. 
    2. Practice implementing a Nextflow process from a bash script.  
    3. Practice chaining Nextflow processes with channels.  
    4. Understand the components of a process such as `input`, `output`,
    `script`, `directive`, and the `workflow` scope.  

In this section you will implement the `QUANTIFICATION` process. This section
contains more self-directed exercises for you to practice applying the
concepts introduced in the previous sections.

![](img/3.excalidraw.png)

You will implement the `02_quant.sh` bash script as a Nexflow process called
`QUANTIFICATION`.  

```bash title="02_quant.sh"
SAMPLE_ID=gut
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"

salmon quant --libType=U -i results/salmon_index -1 ${READS_1} -2 ${READS_2} -o results/${SAMPLE_ID}
```

- `--libType=U` is a required argument. Leave this as is for the script definition.  
- `-i results/salmon_index` is the directory output by the `INDEX` process.  
- `-1` and `-2` are flags for the respective paired reads (`.fq`).  
- `-o` indicates that the command will output files into a directory called `results/gut`

!!! Tip

    The following exercises break down each step in the same order we have been
    following when implementing a process: 

    1. Start with the empty `process` scaffold.  
    2. Define the `script` based on the bash script `02_quant.sh`.  
    3. Identify and define the `output`(/s).  
    4. Identify and define the `input`(/s).  
    5. Define the `directive`s (`container`, `publishDir`).  
    6. Add the `workflow` scope and assign to a channel.  
    7. Create channels and input to the process.  

## 2.3.1 Implementing the process  

### 1. The empty process  

Here is the empty `process` template to get you started:  

```groovy title="main.nf"
process QUANTIFICATION {
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

### 2. Define the process `script`  

!!! question "Exercise"

    Add the `script` definition. Remember to convert hardcoded paths and bash
    variables to Nextflow `process` variables.  

    ??? note "Solution"

        ```groovy title="main.nf"
        process QUANTIFICATION {
          [ directives ]
        
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

### 3. Define the process `output`

!!! question "Exercise"

    Define the `output` definition.  

    ??? tip "Hint"  

        The `output` is a directory of `$sample_id`.

    ??? note "Solution"

        ```groovy title="main.nf"
        process QUANTIFICATION {
          [ directives ]
        
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

```groovy title="main.nf"
process QUANTIFICATION {
  [ directives ]

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

```groovy title="main.nf"
process QUANTIFICATION {
  [ directives ]

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

### 5. Add the process `directives`

!!! question "Exercise"

    Add the `container` and `publishDir` directives to the `QUANTIFICATION`
    process. The container is the same as the `INDEX` process.  

    ??? note "Solution"

        ```groovy title="main.nf"
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

### 6. Call the process in the `workflow` scope  

!!! question "Exercise"

    Add `QUANTIFICATION` to the workflow scope, add the inputs, and assign it
    to a channel called `quant_ch`.  

    ??? note "Solution"

        ```groovy title="main.nf"
                fastq_ch = INDEX(params.transcriptome_file)
                quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
        }
        ```

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

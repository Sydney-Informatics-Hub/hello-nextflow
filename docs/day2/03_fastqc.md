# 2.3 Quality control

In this step, you will step through exercises to implement a quality control
process `FASTQC`, and draw on the concepts introduced in the previous sections.

## 2.3.1 Building the process  

!!! question "Exercise"

    Add a workflow scope that:  
        - Runs process `FASTQC`  
        - Takes `read_pairs_ch` as input  
        - Output is assigned to a channel called `fastqc_ch`  

    ??? note "Solution"

        ```groovy title="main.nf"
            quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
            fastqc_ch = FASTQC(read_pairs_ch)
        }
        ```

The command for the script definition will be:  

```bash
fastqc --outdir [outdir] --format fastq [fastq_1] [fastq_2]
```

Convert to a process script definition:  

```groovy title="main.nf"
process FASTQC {

    script:
    """
    fastqc --outdir "fastqc_logs" --format fastq $reads_1 $reads_2
    """
}
```

!!! question "Exercise"

    Add the `input` definition.

    ??? note "Solution"

        ```groovy title="main.nf"
        process FASTQC {
       
            input:
            tuple val(sample_id), path(reads_1), path(reads_2)
            
            script:
            """
            fastqc --outdir "fastqc_logs" --format fastq $reads_1 $reads_2
            """
        }
        ```

!!! question "Exercise"

    Add the `output` definition.

    ??? note "Solution"

        ```groovy title="main.nf"
        process FASTQC {
       
            input:
            tuple val(sample_id), path(reads_1), path(reads_2)
            
            output:
            path "fastqc_logs"

            script:
            """
            fastqc --outdir "fastqc_logs" --format fastq $reads_1 $reads_2
            """
        }
        ```

!!! question "Exercise"

    Add the `container` and `publishDir` directives. The container for
    `fastqc` is `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0`.

    ??? note "Solution"

        ```groovy title="main.nf"
        process FASTQC {
       
            container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
            publishDir params.outdir, mode: 'copy'

            input:
            tuple val(sample_id), path(reads_1), path(reads_2)
            
            output:
            path "fastqc_logs"

            script:
            """
            fastqc --outdir "fastqc_logs" --format fastq $reads_1 $reads_2
            """
        }
        ```

> Some sort of progress checkpoint here would be good

## 2.3.2 Demarcating the process with `sample_id`  

> Prose  

!!! question "Exercise"

    Update the relevant process definitions so that `fastqc_gut_logs` is
    output.

    ??? note "Solution"

        ```groovy title="main.nf"
        process FASTQC {
       
            container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
            publishDir params.outdir, mode: 'copy'

            input:
            tuple val(sample_id), path(reads_1), path(reads_2)
            
            output:
            path "fastqc_${sample_id}_logs"

            script:
            """
            fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2
            """
        }
        ```

!!! question "Exercise"

    Add a `tag` directive that indicates the name of the tool and sample_id.  

    ??? note "Solution"

        ```groovy title="main.nf"
        process FASTQC {
       
            input:
            tuple val(sample_id), path(reads_1), path(reads_2)
            
            output:
            path "fastqc_${sample_id}_logs"

            script:
            """
            fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2
            """
        }
        ```

Run the workflow:  

```bash
nextflow run main.nf -resume
```

!!! question "Poll"

    What was the command error recevied? Ideas to resolve?

    ??? note "Solution"

        The error:  

        ```console title=main.nf
        Specified output directory 'fastqc_gut_logs' does not exist
        ```

        Fix by updating the script definition:  

        ```groovy title="main.nf"
        process FASTQC {
       
            input:
            tuple val(sample_id), path(reads_1), path(reads_2)
            
            output:
            path "fastqc_${sample_id}_logs"

            script:
            """
            mkdir "fastqc_${sample_id}_logs"
            fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2
            """
        }
        ```

> Inspect .html and ask Poll about it?  


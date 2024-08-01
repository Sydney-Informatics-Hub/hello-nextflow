# hello-nextflow
Training materials for a Nextflow beginners workshop 2024

## Day 2  

Adapated from
[nf-training](https://github.com/nextflow-io/training/blob/master/nf-training/script7.nf).  

### Installation (temporary)  
```bash
mamba create -n day2
mamba activate day2
mamba install -c bioconda nextflow salmon # fastqc, multiqc in docker
```

### Usage  
```bash
nextflow run main.nf
```

Finished running within seconds on laptop with specs:
> CPU: 12th Gen Intel i7-1265U (12) @ 4.800GHz  
> Memory: ~32GB

Note: `*.nf` symlinks to `scripts/*.nf` to simulate running `main.nf` in `projectDir`.  


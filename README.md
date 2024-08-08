# hello-nextflow
Training materials for a Nextflow beginners workshop 2024

## Day 2  

Adapated from
[nf-training](https://github.com/nextflow-io/training/blob/master/nf-training/script7.nf).  

### Installation (temporary)  

Install miniconda: 
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
~/miniconda3/bin/conda init
source ~/.bashrc
```

Install mamba: 
```bash
conda create -n mamba-env -c conda-forge mamba
conda activate mamba-env
mamba --version #to confirm
```

Activate the environment:
```bash
mamba create -n day2
mamba activate day2
mamba install -c bioconda nextflow salmon # fastqc, multiqc in docker  

# docker

# mkdocs
mamba install -c conda-forge mkdocs mkdocs-material
```

### Usage  
```bash
nextflow run main.nf
```

Finished running within seconds on laptop with specs:
> CPU: 12th Gen Intel i7-1265U (12) @ 4.800GHz  
> Memory: ~32GB

**Note:** Numbered `.nf` files are for checkpoints throughout the workshop
when the workflow needs to be run. Participants will build on the one
`main.nf`.  

### mkdocs  

```bash
# mkdocs new .
mkdocs build
```

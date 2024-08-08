# hello-nextflow
Training materials for a Nextflow beginners workshop 2024

## Day 2  

Adapated from
[nf-training](https://github.com/nextflow-io/training/blob/master/nf-training/script7.nf).  

### Installation (temporary)  

Install mamba and packages:
```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
export PATH="$HOME/miniforge3/bin:$PATH"
source ~/.bashrc
conda install mamba -n base -c conda-forge
mamba --version #to confirm

mamba create -n day2
mamba activate day2
mamba install -c bioconda nextflow salmon # fastqc, multiqc in docker  

# mkdocs
mamba install -c conda-forge mkdocs mkdocs-material
```

Install docker:
```bash
# docker https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

sudo docker run hello-world

sudo chmod 666 /var/run/docker.sock
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

To generate pretty html docs:
```bash
cd ~/hello-nextflow/
mkdocs serve
```

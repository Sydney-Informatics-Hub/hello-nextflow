# hello-nextflow
Training materials for a Nextflow beginners workshop 2024

## Day 2  

Adapated from
[nf-training](https://github.com/nextflow-io/training/blob/master/nf-training/script7.nf).

### Installation (temporary)  

#### `mamba`  

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
# complete prompts with defaults

# confirm installation versions
mamba --version
```

```console
mamba 1.5.8
conda 24.3.0
```

Install packages:  

```bash
mamba create -n day2
mamba activate day2
mamba -c bioconda nextflow # process tools via docker containers
mamba -c conda-forge mkdocs mkdocs-material
```

```bash
mamba list | grep -E "nextflow|mkdocs |mkdocs-material"
```

```console
mkdocs                    1.6.0              pyhd8ed1ab_0    conda-forge
mkdocs-material           9.5.31             pyhd8ed1ab_0    conda-forge
mkdocs-material-extensions 1.3.1              pyhd8ed1ab_0    conda-forge
nextflow                  24.04.4              hdfd78af_0    bioconda
```

#### Docker   

Follows [ubuntu installation](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository) and [linux post-install steps](https://docs.docker.com/engine/install/linux-postinstall/).

```bash
# Add Docker's official GPG key:
sudo apt update
sudo apt install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt update

# Install the latgest Docker packages
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Check installed correctly
sudo docker run hello-world

# post-install mods for non-root permissions
sudo group add docker
sudo usermod -aG docker $USER
newgrp docker

# confirm user access
docker run hello-world
```

##### Pull containers  

```bash
docker pull quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0
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

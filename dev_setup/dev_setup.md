
# Developer Installation  

## `mamba`  

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
# complete prompts with defaults

export PATH="$HOME/miniforge3/bin:$PATH"
source ~/.bashrc

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
mamba install -c bioconda nextflow # process tools via docker containers
mamba install -c conda-forge mkdocs mkdocs-material
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

## Docker   

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
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker

# confirm user access
docker run hello-world
```

### Pull containers  

```bash
docker pull quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0
docker pull quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
docker pull quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0
```

## Developer Usage

### VM testing  

#### Part 1
```
nextflow run nextflow-io/hello
```

#### Part 2

**Note:** `.main.nf` and `.nextflow.config` runs the final pipeline. 

```bash
git clone https://github.com/Sydney-Informatics-Hub/hello-nextflow.git
cd hello-nextflow/part2  
mv .nextflow.config nextflow.config
```

Run completed pipeline (includes introspection reports etc., multi-sample and multiple cpus):  

```bash
nextflow run .main.nf --reads data/samplesheet_full.csv
```

## mkdocs  

To generate html docs during development (before deploying):

```bash
cd ~/hello-nextflow/
mkdocs serve
# open http://127.0.0.1:8000/ in browser
```

To render docs for website (manually):  

```bash
mkdocs gh-deploy
```

Open [https://sydney-informatics-hub.github.io/hello-nextflow/](https://sydney-informatics-hub.github.io/hello-nextflow/) in a browser.  


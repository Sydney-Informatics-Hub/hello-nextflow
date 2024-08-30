# hello-nextflow
Training materials for a Nextflow beginners workshop 2024

## Developer Installation  

### `mamba`  

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

### Docker

Follow Docker instructions under "User" Installation.  

### mkdocs  

To render docs for website:  

```bash
# mkdocs new .
mkdocs build
```

To generate html docs during development:
```bash
cd ~/hello-nextflow/
mkdocs serve
# open http://127.0.0.1:8000/ in browser
```

## "User" Installation  

Installing dependencies on account with root access for all (non-root) users.  

### Java  

Use `sdkman` for java install, recommended by 
[Nextflow docs](https://www.nextflow.io/docs/latest/install.html).

```bash
# sdkman dependencies
cd $HOME
sudo apt update
sudo apt upgrade -yU
sudo apt install zip unzip tree # tree for day1

# Install latest java LTS version of Temurin with sdkman
curl -s https://get.sdkman.io | bash
source ~/.bashrc
sdk install java 17.0.10-tem

# Move for access for all users  
sudo mkdir /usr/lib/jvm
sudo mv ${SDKMAN_CANDIDATES_DIR}/java/17.0.10-tem /usr/lib/jvm
sudo ln -s /usr/lib/jvm/17.0.10-tem/bin/java /usr/bin/java
```

Validate `java` install:  

```bash
java -version
```

```console
openjdk version "17.0.10" 2024-01-16
OpenJDK Runtime Environment Temurin-17.0.10+7 (build 17.0.10+7)
OpenJDK 64-Bit Server VM Temurin-17.0.10+7 (build 17.0.10+7, mixed mode, sharing)
```

### Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod 755 nextflow
sudo mv nextflow

# Validate
nextflow info
```

```console
  Version: 24.04.4 build 5917
  Created: 01-08-2024 07:05 UTC
  System: Linux 6.8.0-35-generic
  Runtime: Groovy 4.0.21 on OpenJDK 64-Bit Server VM 21-internal-adhoc.conda.src
  Encoding: UTF-8 (UTF-8)
```

### Docker   

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

#### Pull containers  

```bash
docker pull quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0
docker pull quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
docker pull quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0
```

## Usage (testing)   

```bash
git clone https://github.com/Sydney-Informatics-Hub/hello-nextflow.git # need some extra args for draft branch  
cd hello-nextflow/day2  
```

Complete metadata and config to run the full, final pipeline:  

```bash
echo "docker.enabled=true" > nextflow.config
echo "liver,data/ggal/liver_1.fq,data/ggal/liver_2.fq" >> data/samplesheet.csv
echo "lung,data/ggal/lung_1.fq,data/ggal/lung_2.fq" >> data/samplesheet.csv
```

Run:  

```bash
nextflow run main.nf
```

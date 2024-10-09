# VM set up 

## VM Installation  

## Install script  

**Things to include:**  

- Logging  
- Setting directory permissions  
- Checking required network access  
- Check available disk  
- Automated testing to ensure things are set up correctly, e.g.  
    - Permissions  
    - Install success  
    - Use 0 vs non-0 exit status  
- Automated cleanup  
- Setting and verifying environmental variables  
- Checking data integrity (checksums)  

## Usage  

```bash
cd $HOME
git clone https://github.com/Sydney-Informatics-Hub/hello-nextflow.git

# Sort out docker groups before installation 
uname=user1
sudo adduser $uname
if ! getent group docker; then
    sudo groupadd docker
fi
sudo usermod -aG $USER
sudo usermod -aG $uname
newgrp docker # not sure if this is needed if rebooting
sudo reboot

# Install and validate for root user, set docker permissions for user
bash hello-nextflow/vm_install/ubuntu_install.sh $uname

# Test installs for user
sudo su - $uname

## utils
tree
tr

## nextflow
nextflow run nextflow-io/hello # in $HOME
rm -rfv .nextflow* work/

## docker
docker run "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
docker run "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
docker run "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"

# Set dir structure and get files required for the workshop (in user1)
cd $HOME
git clone https://github.com/Sydney-Informatics-Hub/hello-nextflow.git
cp -r hello-nextflow/part1 hello-nextflow/part2
rm -rfv hello-nextflow/
```

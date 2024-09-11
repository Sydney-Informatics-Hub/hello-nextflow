#!/usr/bin/bash -eu

uname=$1 # For extra (non-root) user

# Set timezone (for logging)  
TZ=$(timedatectl show --property=Timezone --value)
if [ $TZ != "Australia/Sydney" ]; then
    echo "Setting timezone to Australia/Sydney"
    sudo timedatectl set-timezone Australia/Sydney
fi

# Define logger
log() {
    TIMESTAMP=`date +%y-%m-%d\ %H:%M:%S`
    echo "[${TIMESTAMP}] $@"
}

# Check executable
check_executable() {
    if [ ! -x "$(command -v $1)" ]; then
        log "ERROR: $1 not executable"
        exit 1
    fi
}

## apt
log "Updating and upgrading packages..."
sudo apt update -y && sudo apt upgrade -y

log "Installing required system utilities..."
SYSUTILS=(zip unzip tree curl)
sudo apt install -y "${SYSUTILS[@]}"

log "Checking system utilities are executable..."
for util in "${SYSUTILS[@]}"; do
    check_executable $util
done
check_executable tr

log "Installing default-jre"
sudo apt install default-jre

log "Verifying java install..."
check_executable "java -version"

## Nextflow
log "Installing nextflow..."
curl -s https://get.nextflow.io | bash
chmod 755 ./nextflow
sudo mv ./nextflow /usr/local/bin

log "Verifying nextflow install..."
check_executable nextflow

## Docker   
# Follows ubuntu installation docs and linux post-install steps
# https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository 

log "Adding Docker's official GPG key..."
sudo apt update -y
sudo apt install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

log "Adding the Docker repository to Apt sources..."
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt update -y

log "Installing the latest Docker packages..."
sudo apt install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

log "Verifying docker install for $USER (root)..."
if ! sudo docker run hello-world; then
    log "ERROR: Failed to run `sudo docker run hello-world`"  
    exit 1
fi

log "Verifying docker install for $USER (non-root)..."
if ! docker run hello-world; then
    log "ERROR: Failed to run `sudo docker run hello-world`"  
    exit 1
fi

log "Pulling docker containers"
IMAGES=(
	"quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
	"quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
	"quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
)

for image in "${IMAGES[@]}"; do
    docker pull $image
done

log "Validating docker containers"
for image in "${IMAGES[@]}"; do
    if ! docker run $image; then
       	log "ERROR: Failed to run $image"
    	exit 1
    fi
done

log "Installation for $USER successful!"

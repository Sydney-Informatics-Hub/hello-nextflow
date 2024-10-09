# Use the Gitpod workspace base image
FROM gitpod/workspace-base

USER root

# Set environment variables
ENV TZ=Australia/Sydney
ENV DEBIAN_FRONTEND=noninteractive

# Set the timezone and install required utilities
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
    apt-get update --quiet && \
    apt-get install --quiet --yes \
        apt-transport-https \
        apt-utils \
        sudo \
        git \
        less \
        wget \
        curl \
        tree \
        graphviz \
        zip \
        unzip \
        default-jre \
        ca-certificates

# User permissions
RUN mkdir -p /workspace/data \
&& chown -R gitpod:gitpod /workspace/data

# Change user to gitpod
USER gitpod

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    chmod 755 ./nextflow && \
    mv ./nextflow /usr/local/bin

# Create .nextflow directory and set proper permissions
RUN mkdir -p /home/gitpod/.nextflow && \
chown -R gitpod:gitpod /home/gitpod/.nextflow

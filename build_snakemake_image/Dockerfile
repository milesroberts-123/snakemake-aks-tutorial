# Start from ubuntu base image, make sure its updated so packages work, then install basic necessities
FROM ubuntu:20.04

RUN apt-get update --no-install-recommends --assume-yes && apt-get upgrade --no-install-recommends --assume-yes

RUN apt-get install --no-install-recommends --assume-yes wget curl bzip2 ca-certificates gnupg2 squashfs-tools git

# Add opt/conda to environment path
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Install mamba to speed up conda
# I needed to remove the strict channel priorities from the system config to get all of the packages installed. I'm not sure if this will cause future problems
RUN curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh > mambaforge.sh
RUN bash mambaforge.sh -b -p /opt/conda
RUN conda config --system
RUN rm mambaforge.sh

# copy conda environment into docker container
ADD ./snakemake-docker-environment.yaml .
RUN mamba env update --file ./snakemake-docker-environment.yaml && conda clean -tipy

# activate snakemake environment upon starting container, so that installed software is accessible to kubernetes
RUN echo "source activate snakemake" > ~/.bashrc
ENV PATH=/opt/conda/envs/snakemake/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

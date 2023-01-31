# Ubuntu base image, make sure its updated so packages work
FROM ubuntu:20.04

RUN apt-get update --no-install-recommends --assume-yes && apt-get upgrade --no-install-recommends --assume-yes

RUN apt-get install --no-install-recommends --assume-yes wget curl bzip2 ca-certificates gnupg2 squashfs-tools git

# Install miniconda
# ENV PATH $CONDA_DIR/bin:$PATH
# RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda$CONDA_PYTHON_VERSION-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
 #   echo 'export PATH=$CONDA_DIR/bin:$PATH' > /etc/profile.d/conda.sh && \
 #   /bin/bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
 #   rm -rf /tmp/*

# install mamba to speed up conda
# RUN conda install -y mamba -c conda-forge

# Add opt/conda to environment path
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Install mamba
RUN curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh > mambaforge.sh
RUN bash mambaforge.sh -b -p /opt/conda
RUN conda config --system
RUN rm mambaforge.sh

# copy conda environment into docker file
# Add install path for where kubernetes seems to be looking
ADD ./environment.yml .
RUN mamba env update --file ./environment.yml && conda clean -tipy

# activate snakemake environment upon start up
RUN echo "source activate snakemake" > ~/.bashrc
ENV PATH=/opt/conda/envs/snakemake/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

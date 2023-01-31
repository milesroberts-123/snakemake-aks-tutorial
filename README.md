# Running snakemake using Azure kubernetes service

This walkthrough works as of 2023-01-31. Unfortunately, future changes to either Azure Kubernetes Service or snakemake may prevent this walkthrough from working. However, this will serve as a record of my debugging efforts to help any poor soul that attempts using Azure Kubernetes. Good luck and may Bill Gates have mercy on you.

## Contents

[Main changes from official snakemake tutorial](#main-changes-from-official-snakemake-tutorial)

[References](#references)

0. [Set up linux environment to lauch kubernetes](#set-up-linux-environment-to-lauch-kubernetes)

1. [Create a conda environment with needed dependencies](#Create-a-conda-environment-with-needed-dependencies)

2. [Create a storage account](#Create-a-storage-account)

3. [Create auto-scaling kubernetes cluster](#create-auto-scaling-kubernetes-cluster) 

4. [Optionally build our own docker image of snakemake](#optionally-build-our-own-docker-image-of-snakemake)

5. [Run the workflow](#run-the-workflow)

## Main changes from official snakemake tutorial

A snakemake tutorial repo can be found [here](https://github.com/snakemake/snakemake-tutorial-data), but it doesn't work for AKS out of the box. 

I modified the above resource in the following ways:

* copied [this workflow](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#summary) to a Snakefile

* Added a `scripts/plot-quals.py` script from [here](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-6-using-custom-scripts)

* Removed `data/samples/C.fastq` because it was unecessary to provide a minimal snakemake example and also not included in Snakefile

* Added `conda: "enviornment.yaml"` to each snakemake rule so that kubernetes would install software on the fly

* Added genome BWA index to input of the `bwa_map` rule. BWA mem needs the index in order to run, but kubernetes won't download the index to its nodes unless the index is defined in the input.

## References

The material in this walkthrough was compiled from the following links:

[azure-specific snakemake+kubernetes tutorial](https://snakemake.readthedocs.io/en/stable/executor_tutorial/azure_aks.html)

[generic snakemake+kubernetes tutorial](https://snakemake.readthedocs.io/en/stable/executing/cloud.html)

[install azure-cli](https://learn.microsoft.com/en-us/cli/azure/install-azure-cli-linux?pivots=script)

[Dockerfile for latest snakemake image](https://hub.docker.com/layers/snakemake/snakemake/latest/images/sha256-6e33aee97fbca79e99a9603ba7502ef732c466bf96b912da9c390189431f04ea?context=explore)

[How to replicate conda environment in a Dockerfile](https://blog.ceshine.net/post/replicate-conda-environment-in-docker/#export-conda-environment)

[Video of how to build image and push it to Docker hub](https://www.youtube.com/watch?v=iqqDU2crIEQ&t=1002s)

And a few github issues:

* https://stackoverflow.com/questions/68681883/snakemake-auto-scaling-azure-kubernetes-cluster-without-shared-filesystem-fail

* https://github.com/snakemake/snakemake/issues/1167

## Set up linux environment to lauch kubernetes

First we need a computer from which to lauch kubernetes. I'm running windows locally and the tutorials above suggest I need the `azure-cli` installed. I'm guessing my windows subsystem for linux would cause problems for this, so I'm thinking I should just get an azure data science virtual machine (which includes azure-cli) in the cloud to save myself the hassel.

### create linux vm in azure portal

Find data science virtaul machine under All services > Create a resource > Marketplace, then use the following specifications:

resource group: ccf22_robe1195 (or whatever you want)

vitual machine name: ccf22robe1195snakemaketest1 (or whatever you want)

image: data science virtual machine ubuntu 20.04

size: b2s, 2 cpus, 4 GB

username: robe1195 (or whatever you want)

password: **************************************** (40 characters, randomly generated)

os disk type: standard ssd

tag name: mainproject

value: test

### log into vm from local terminal (for me: a windows subsystem for linux)

Get IP address by clicking on the VM resource

Log on and enter password with `ssh robe1195@<insert IP address>`. You may need to wait a bit before connection is accepted

Check if conda is installed with `conda --help`

Check if azure-cli is installed with `az --help`

Check if docker is installed with `docker --help`

## Create a conda environment with needed dependencies

Once you log into your VM, you need to get snakemake, azure-blob-storage, and kubernetes modules all in the same environment.

First create the conda environment and install snakemake:

```
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

Now install other dependencies into the snakemake environment

```
# activate snakemake environment
conda activate snakemake

# install kubernetes
# For some reason, installing via both pip and conda prevents later issues
conda install -c conda-forge kubernetes
pip install kubernetes

# install software for interacting with blob storage
conda install -c conda-forge azure-storage-blob
```

Again, Azure CLI should already be installed because I'm using the DSVM, so that saves one huge headache

## Create a storage account

We want snakemake and kubernetes to write to blob storage, so we first need to create a storage account. 

You could do this manually in the portal too or use the following commands

```
# change the following names as required
# azure region where to run:
region=northcentralus
# name of the resource group to create (just use the same one assigned to me during the fellowship)
resgroup=ccf22_robe1195
# name of storage account to create (all lowercase, no hyphens etc.):
stgacct=ccf22robe1195snakekub

# create a resource group with name and in region as defined above
# I don't need to do this because I already have a resource group under the cloud fellowship subscription
# az group create --name $resgroup --location $region

# loging to azure, run the below code then follow the steps
az login --use-device-code

# create a general purpose storage account with cheapest SKU
az storage account create -n $stgacct -g $resgroup --sku Standard_LRS -l $region

# get storage account key for later use
# the below line is in the tutorial, but doesn't work and the tutorial contains a typo ($storageacct -> $stgacct)
# stgkey=$(az storage account keys list -g $resgroup -n $stgacct | head -n1 | cut -f 3)
# Instead we'll use json file processing and remove quotes to get the key
stgkey=$(az storage account keys list -g $resgroup -n $stgacct | jq .[1].value | sed 's/^"//' | sed 's/"$//')

# finally, create the storage container
az storage container create --resource-group $resgroup --account-name $stgacct --account-key $stgkey --name snakemake-tutorial

## Get example data from this repo

Download this repo to your VM, then upload the `data/` portion to your storage account so that kubernetes can access it

```
mkdir tutorial
cd tutorial
git https://github.com/milesroberts-123/snakemake-aks-tutorial.git
cd snakemake-aks-tutorial
az storage blob upload-batch -d snakemake-tutorial --account-name $stgacct --account-key $stgkey -s data/ --destination-path data
```

## Create auto-scaling kubernetes cluster

```
# change the cluster name as you like
# needed to add --generate-ssh-keys to tutorial command
clustername=snakemaks-aks
az aks create --generate-ssh-keys --resource-group $resgroup --name $clustername --vm-set-type VirtualMachineScaleSets --load-balancer-sku standard --enable-cluster-autoscaler --node-count 1 --min-count 1 --max-count 3 --node-vm-size Standard_B4ms
```

Now fetch credentials for cluster so that you can interact with it

```
# get credentials
az aks get-credentials --resource-group $resgroup --name $clustername

# print basic cluster info
kubectl cluster-info
```

## Optionally build our own docker image of snakemake

I already did this step and pushed the result to dockerhub, so you can technically skip it

An important lesson I learned is that the snakemake version on which your lauch vm needs to match the snakemake container used in kubernetes. However, a problem with the current snakemake docker image is that it excludes the azure-blob-storage dependency. Thus, I built a new snakemake image that has this dependency installed.

I started by exporting my current snakemake environment, which already had the needed dependencies (snakemake, azure-blob-storage, and kubernetes)

```
cd build_snakemake_image
conda env export > snakemake-docker-environment.yaml
```

Then copied the environment over to a Dockerfile that looks like this:

```
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
```

Now build the docker image

`docker build -t ccf22robe1195snakemake:latest .`

Now I needed to push the image to dockerhub after logging into docker

First I needed to make the `snakemake-aks` repo on docker hub, then execute:

```
docker tag ccf22robe1195snakemake milesroberts/snakemake-aks
docker logout
docker login
docker push milesroberts/snakemake-aks
```

Finally return to the directory with the Snakefile with `cd ..`

## Run the workflow

Export storage account and key information, then call snakemake with `--kubernetes` using our custom snakemake image

```
export AZ_BLOB_ACCOUNT_URL="https://${stgacct}.blob.core.windows.net"

export AZ_BLOB_CREDENTIAL="$stgkey"

snakemake --kubernetes --container-image docker.io/milesroberts/snakemake-aks:latest --default-remote-prefix snakemake-tutorial --default-remote-provider AzBlob --envvars AZ_BLOB_ACCOUNT_URL AZ_BLOB_CREDENTIAL --use-conda --jobs 3
```

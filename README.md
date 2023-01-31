# Running snakemake using Azure kubernetes service

This walkthrough works as of 2023-01-31. Unfortunately, future changes to either Azure Kubernetes Service or snakemake may prevent this walkthrough from working. However, this will serve as a record of my debugging efforts to help any poor soul that attempts using Azure Kubernetes.

## Contents

[Main changes from official snakemake tutorial](#main-changes-from-official-snakemake-tutorial)

[References](#references)

0.[Set up linux environment to lauch kubernetes](#set-up-linux-environment-to-lauch-kubernetes)

## Main changes from official snakemake tutorial

A snakemake tutorial repo can be found [here](https://github.com/snakemake/snakemake-tutorial-data), but it doesn't work for AKS out of the box. 

I modified the above resource in the following ways:

* copied [this workflow](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#summary) to a Snakefile

* Added a `scripts/plot-quals.py` script from [here](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-6-using-custom-scripts)

* Removed `data/samples/C.fastq` because it was unecessary to provide a minimal snakemake example and also not included in Snakefile

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

I'm running windows locally and the tutorials above suggest I need the `azure-cli` installed. I'm guessing my windows subsystem for linux would cause problems for this, so I'm thinking I should just get an azure data science virtual machine (which includes azure-cli) to save myself the hassel.

### create linux vm in azure portal

Find data science virtaul machine under All services > Create a resource > Marketplace, then use the following specifications

resource group: ccf22_robe1195
vitual machine name: ccf22robe1195snakemaketest1
image: data science virtual machine ubuntu 20.04
size: gonna use a super cheap standard b2s, 2 cpus, 4 GB
username: robe1195
password: **************************************** (40 characters, randomly generated)
os disk type: standard ssd
tag name: mainproject
value: test2

### log into vm from local terminal (for me: a windows subsystem for linux)

Get IP address by clicking on the VM resource

Log on and enter password with `ssh robe1195@<insert IP address>`. May need to wait a bit before connection is accepted

Check if conda is installed with `conda --help`

Check if azure-cli is installed with `az --help`

Check if docker is installed with `docker --help`

## 2. Create a conda environment

To create the conda environment:

```
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

Now install other dependencies into the snakemake environment
```
conda activate snakemake

# install kubernetes
conda install -c conda-forge kubernetes

pip install kubernetes

# install software for interacting with blob storage
conda install -c conda-forge azure-storage-blob
```

Azure CLI should already be installed because I'm using the DSVM, so that saves one huge headache

## 3. Create a storage account

I could do this manually in the portal too, but I think I'll try the commands to see how it goes. I need to learn how to automate things eventually. Might as well start now!

```
# change the following names as required
# azure region where to run:
region=northcentralus
# name of the resource group to create (just use the same one assigned to me during the fellowship)
resgroup=ccf22_robe1195
# name of storage account to create (all lowercase, no hyphens etc.):
stgacct=ccf22robe1195snakekub2

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
# instead just run
az storage account keys list -g $resgroup -n $stgacct

# then copy and paste one of the values to an environment variable
# stgkey=<value printed from above command>
# or use json file processing and remove quotes
stgkey=$(az storage account keys list -g $resgroup -n $stgacct | jq .[1].value | sed 's/^"//' | sed 's/"$//')

# create the storage container
az storage container create --resource-group $resgroup --account-name $stgacct --account-key $stgkey --name snakemake-tutorial

## 4. Upload data to storage account

```
mkdir tutorial
cd tutorial
git clone https://github.com/snakemake/snakemake-tutorial-data.git
cd snakemake-tutorial-data
az storage blob upload-batch -d snakemake-tutorial --account-name $stgacct --account-key $stgkey -s data/ --destination-path data
```

## 5. Create auto-scaling kubernetes cluster

First make the cluster

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

## 6. get snakefile

https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#summary

Copy this workflow into a snakefile in the snakemake-tutorial directory. To the samples at the top, add "C" if desired

Create plot-quals.py script from here: https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-6-using-custom-scripts

```
mkdir scripts/
cd scripts/
nano plot-quals.py
```

## try building our own docker image of snakemake, include azure software, from conda environment

export conda environment

`conda env export > environment.yml`

write this dockerfile, using this blog as a guide: https://blog.ceshine.net/post/replicate-conda-environment-in-docker/#export-conda-environment

```
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
ADD ./environment.yml .
RUN mamba env update --file ./environment.yml && conda clean -tipy
```

Now build image

`docker build -t ccf22robe1195snakemake:latest .`

Now push the image to dockerhub after logging into docker

First needed to make snakemake-aks repo on docker hub

This video was helpful for this part: https://www.youtube.com/watch?v=iqqDU2crIEQ&t=1002s

```
docker tag ccf22robe1195snakemake milesroberts/snakemake-aks
docker logout
docker login
docker push milesroberts/snakemake-aks
```

## 6. run the workflow


```
export AZ_BLOB_ACCOUNT_URL="https://${stgacct}.blob.core.windows.net"

export AZ_BLOB_CREDENTIAL="$stgkey"

snakemake --kubernetes --container-image docker.io/snakemake/snakemake:latest --default-remote-prefix snakemake-tutorial --default-remote-provider AzBlob --envvars AZ_BLOB_ACCOUNT_URL AZ_BLOB_CREDENTIAL --use-conda --jobs 3


I seem to get the following issue:

```
$ kubectl logs snakejob-cf7f1188-70e9-5f17-8a96-8d6a65b62fa7
WorkflowError:
The Python 3 package 'azure-storage-blob' need to be installed to use Azure Storage remote() file functionality. No module named 'azure'
  File "/opt/conda/envs/snakemake/lib/python3.11/importlib/__init__.py", line 126, in import_module
  File "<frozen importlib._bootstrap>", line 1206, in _gcd_import
  File "<frozen importlib._bootstrap>", line 1178, in _find_and_load
  File "<frozen importlib._bootstrap>", line 1149, in _find_and_load_unlocked
  File "<frozen importlib._bootstrap>", line 690, in _load_unlocked
  File "<frozen importlib._bootstrap_external>", line 940, in exec_module
  File "<frozen importlib._bootstrap>", line 241, in _call_with_frames_removed
```

I think this is because the azure module was deprecated and now snakemake should be using az or something similar

If I use the docker.io/snakemake/snakemake:v6.1.1 container image, I get the error:

`azure.core.exceptions.ResourceNotFoundError: The specified container does not exist.`

I guess now the key is to find the correct snakemake container

Or maybe I need to modify the container image url. The default value for `--container-image` is https://hub.docker.com/r/snakemake/snakemake

### after building my own docker image of snakemake, I tried this:

snakemake --kubernetes --container-image docker.io/milesroberts/snakemake-aks:latest --default-remote-prefix snakemake-tutorial --default-remote-provider AzBlob --envvars AZ_BLOB_ACCOUNT_URL AZ_BLOB_CREDENTIAL --use-conda --jobs 3

but then I got this error

`/opt/conda/bin/python: No module named snakemake`

so I changed the prefix in my environment.yml file to `prefix: /opt/conda`, rebuilt the docker image, and re-pushed it to docker hub

I then ran the same snakemake command as above using my own docker image `milesroberts/snakemake-aks` and got the same error

Now I'm gonna try specifying the prefix in the last command of the dockerfile

`RUN mamba env update --file ./environment.yml -p /opt/conda && conda clean -tipy`

but then got these two really obscure errors after everything installed successfully

`OSError: Could not find a suitable TLS CA certificate bundle, invalid path: /opt/conda/lib/python3.10/site-packages/certifi/cacert.pem`

`ModuleNotFoundError: No module named 'conda.cli.main_info'`

So then I just specified for snakemake environment to activate upon boot in the bashrc file

rebuild, retag, repush the resulting image

snakemake then appears to run! but I got this error

```
/usr/bin/bash: bwa: command not found
/usr/bin/bash: samtools: command not found

```

So then to the snakemake rules I specified the environment.yaml file that came with the snakemake tutorial with

`conda: environment.yaml`

This worked! However, bwa mem couldn't see the index files in blob storage. I think you need to specify the index files in the snakemake rule in order for them to be downloaded to the job node. Thus, I added the index to bwa_map rule

```
rule bwa_map:
    input:
        genome="data/genome.fa",
        sample="data/samples/{sample}.fastq",
        index1="data/genome.fa.amb",
        index2="data/genome.fa.ann",
        index3="data/genome.fa.bwt",
        index4="data/genome.fa.fai",
        index5="data/genome.fa.pac",
        index6="data/genome.fa.sa"
    output:
        "mapped_reads/{sample}.bam"
    conda: "environment.yaml"
    shell:
        "bwa mem {input.genome} {input.sample} | samtools view -Sb - > {output}"
```

THIS FREAKING WORKS!!!!!!!!!!!!!!

I think my next step is to create an updated snakemake tutorial repo with a snakefile, a plot-quals.py script, and an updated snakemake file

FROM continuumio/miniconda3

# File Author / Maintainer
MAINTAINER Zhixia XIAO <zhixia.xiao@link.cuhk.edu.hk>

## Set root user to install things
USER root
## Set WORKDIR for create conda env
WORKDIR /app

## Change default shell to bash
SHELL ["/bin/bash", "--login", "-c"]

COPY miRDP2_conda.yml .
RUN conda env create -f miRDP2_conda.yml

## Copy rfam_ncRNA index to the image
COPY index/rfam_index.* /opt/conda/envs/miRDP2/bin/scripts/index/

COPY miRDP2-v1.1.4_pipeline.bash entrypoint.sh ./
RUN chmod +x entrypoint.sh
RUN chmod +x miRDP2-v1.1.4_pipeline.bash

## Create user 1000:1000
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 ubuntu
## Set WORKDIR for data analysis, default to /data
WORKDIR /data
## Change permission of /data for user "ubuntu"
RUN chown ubuntu: /data
## Exec the scripts with user ubuntu (uid 1000)
USER ubuntu

## Set entrypoint, with Exec form
ENTRYPOINT ["/app/entrypoint.sh"]
## Default parameter is help
CMD ["-h"]

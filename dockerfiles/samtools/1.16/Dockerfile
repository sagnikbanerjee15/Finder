#######################################################################################################################################################################
# samtools
#######################################################################################################################################################################

FROM condaforge/mambaforge
LABEL maintainer="sagnikbanerjee15@gmail.com" 
LABEL org.opencontainers.image.source https://github.com/sagnikbanerjee15/dockerized_tools_and_pipelines


ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

RUN mamba install -y -c bioconda samtools=1.16.1

#Esyntax=docker/dockerfile:1
FROM ubuntu:20.04
RUN apt update 
RUN DEBIAN_FRONTEND=noninteractive apt install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y r-base file python3-venv libssl-dev libcurl4-openssl-dev python3-pip
RUN R -q -e 'install.packages("BiocManager")' 
RUN R -q -e 'BiocManager::install(c("Seurat", "dplyr", "patchwork", "optparse"))' 
RUN pip3 install umap-learn
WORKDIR /app
COPY . . 
CMD Rscript /app/simple_scrnaseq.R -d /data -o /output

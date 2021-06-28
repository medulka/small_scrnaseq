## scrnaseq

A docker image for running the single cell RNA pipeline Rscript simple_scrnaseq.R. 
The pipeline is inspired by (https://satijalab.org/seurat/) and suitable for analysing droplet-based Genomics 10X data. 
The docker image should work on desirable dataset. You can test it on test_data.tar.gz (added to this repository).

The docker image can be run by:


```bash
docker run -it -d [data] -o [output directory] scrnaseq
```
Please specify 
* your input data directory
* and the desirable output file.


**NOTE:** The docker image can be downloaded from [Docker Hub](https://hub.docker.com/repository/docker/medulka/scrnaseq).



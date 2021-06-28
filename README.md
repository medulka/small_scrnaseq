## scrnaseq_pipeline
simple scrnaseq pipeline based on Seurat

In this repository, you can find a docker image for running the single cell RNA pipeline Rscript simple_scrnaseq.R. 
The pipeline is inspired by https://satijalab.org/seurat/ and suitable for analysing droplet-based Genomics 10X data. 
The docker image should work on desirable dataset. Anyway, you can test it on /test_data/ saved in this repository.

To run pipeline in the Docker, please specify 
* the your data directory
* and the output file where the result will be stored.


```bash
docker run -it -d [data] -o [output directory] scrnaseq
```

You can find the Docker Image on Docker Hub. 

**IMPORTANT:** 

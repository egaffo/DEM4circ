#!/bin/bash

ANN_DIR=/home/enrico/annotation
IDX_DIR=/home/enrico/indexes
#READS_DIR=/home/enrico/datasets/PRJCA000751
READS_DIR=/home/enrico/datasets/ji_2019

# docker run -u `id -u`:`id -g` --rm -it -v $(pwd):/data -v $ANN_DIR:/annotation -v $READS_DIR:/reads -v $IDX_DIR:/indexes egaffo/circompara2:v0.1.2.1 -n
docker run -u `id -u`:`id -g` --rm -it -v $(pwd):/data -v $ANN_DIR:/annotation -v $READS_DIR:/reads -v $IDX_DIR:/indexes egaffo/circompara2:v0.1.2.1 '-j4 -k' 2> ccp.err | tee -a ccp.log

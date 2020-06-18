#!/bin/bash

cd "$(dirname "$0")"

HIGLASS_SERVER_DIR=$1 # move into to the higlass server directory

python $HIGLASS_SERVER_DIR/manage.py delete_tileset \
    --uuid chromsizes-hg19
python $HIGLASS_SERVER_DIR/manage.py ingest_tileset \
    --filename chromSizes.tsv \
    --filetype chromsizes-tsv \
    --datatype chromsizes \
    --coordSystem hg19 \
    --uid chromsizes-hg19

for i in {1..10}
    do
        python $HIGLASS_SERVER_DIR/manage.py delete_tileset \
            --uuid pbmc_10x_peaks_${i}
        python $HIGLASS_SERVER_DIR/manage.py ingest_tileset \
            --filename ./data/processed/pbmc_10x_peaks_${i}.bw \
            --uid pbmc_10x_peaks_${i} \
            --filetype bigwig \
            --datatype vector \
            --coordSystem hg19
done
#!/bin/bash

cd "$(dirname "$0")"
SCRIPT_DIR="$(pwd)"

cd $1 # move into to the higlass server directory

for i in {1..10}
    do
        python manage.py delete_tileset \
            --uuid pbmc_10x_peaks_${i}
        python manage.py ingest_tileset \
            --filename ${SCRIPT_DIR}/data/processed/pbmc_10x_peaks_${i}.bw \
            --uid pbmc_10x_peaks_${i} \
            --filetype bigwig \
            --datatype vector \
            --coordSystem hg19
done
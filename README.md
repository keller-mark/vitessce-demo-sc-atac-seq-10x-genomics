Set up and activate conda environment:

```sh
conda env create -f environment.yml
source activate clodius-env
```

Run snakemake:

```sh
snakemake --cores 1
```

Ingest using higlass-server:

```sh
python manage.py ingest_tileset \
    --filename data/processed/pbmc_10x_peaks_1.bw \
    --filetype bigwig \
    --datatype vector \
    --coordSystem hg19
# for clusters 1 to 10
```
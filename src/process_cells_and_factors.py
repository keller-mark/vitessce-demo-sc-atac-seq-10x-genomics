import json
import pandas as pd

# column names
BARCODE = "Barcode"
TSNE_1 = "TSNE-1"
TSNE_2 = "TSNE-2"
CLUSTER = "Cluster"

if __name__ == "__main__":
    tsne_df = pd.read_csv(snakemake.input["tsne"])
    kmeans_df = pd.read_csv(snakemake.input["kmeans"], index_col=0)

    cells = {}
    factors = {
        "kmeans": {
            "map": [ f"Cluster {i}" for i in range(1, 11) ],
            "cells": {}
        }
    }

    for i, row in tsne_df.iterrows():
        cluster_i = int(kmeans_df.at[row[BARCODE], CLUSTER])
        cells[row[BARCODE]] = {
            "mappings": {
                "t-SNE": [ row[TSNE_1], row[TSNE_2] ]
            },
            "genes": { },
            "xy": [0, 0],
            "factors": {
                "kmeans": f"Cluster {cluster_i}"
            },
            "poly": []
        }
        factors["kmeans"]["cells"][row[BARCODE]] = (cluster_i - 1)
    
    

    with open(snakemake.output["cells"], "w") as f:
        json.dump(cells, f)
    
    with open(snakemake.output["factors"], "w") as f:
        json.dump(factors, f)
import json
import pandas as pd

# column names
BARCODE = "Barcode"
TSNE_1 = "TSNE-1"
TSNE_2 = "TSNE-2"
CLUSTER = "Cluster"

if __name__ == "__main__":
    tsne_1_df = pd.read_csv(snakemake.input["tsne_1"])
    tsne_2_df = pd.read_csv(snakemake.input["tsne_2"], index_col=0)
    kmeans_df = pd.read_csv(snakemake.input["kmeans"], index_col=0)

    cells = {}
    factors = {
        "kmeans": {
            "map": [ f"Cluster {i}" for i in range(1, 11) ],
            "cells": {}
        }
    }

    for i, row in tsne_1_df.iterrows():
        cluster_i = int(kmeans_df.at[row[BARCODE], CLUSTER])
        
        cells[row[BARCODE]] = {
            "mappings": {
                "t-SNE": [ row[TSNE_1], row[TSNE_2] ],
                "t-SNE from LSA": [ tsne_2_df.at[row[BARCODE], TSNE_1], tsne_2_df.at[row[BARCODE], TSNE_2] ]
            },
            "genes": { },
            "xy": [0, 0],
            "factors": {
                "kmeans": f"Cluster {cluster_i}"
            },
            "poly": []
        }
        factors["kmeans"]["cells"][row[BARCODE]] = (cluster_i - 1)
    
    cell_sets = {
        "datasetId": "atac-10x-pbmc-2019",
        "setsType": "cell",
        "version": "0.0.24",
        "setsTree": [
            {
                "key": "kmeans",
                "name": "k-means with 10 clusters",
                "color": [128, 128, 128]
            },
            {
                "key": "kmeans\tCluster 1",
                "name": "Cluster 1",
                "color": [166, 206, 227]
            },
            {
                "key": "kmeans\tCluster 2",
                "name": "Cluster 2",
                "color": [227, 26, 28]
            },
            {
                "key": "kmeans\tCluster 3",
                "name": "Cluster 3",
                "color": [51, 160, 44]
            },
            {
                "key": "kmeans\tCluster 4",
                "name": "Cluster 4",
                "color": [251, 154, 153]
            },
            {
                "key": "kmeans\tCluster 5",
                "name": "Cluster 5",
                "color": [253, 191, 111]
            },
            {
                "key": "kmeans\tCluster 6",
                "name": "Cluster 6",
                "color": [31, 120, 180]
            },
            {
                "key": "kmeans\tCluster 7",
                "name": "Cluster 7",
                "color": [178, 223, 138]
            },
            {
                "key": "kmeans\tCluster 8",
                "name": "Cluster 8",
                "color": [240, 147, 146]
            },
            {
                "key": "kmeans\tCluster 9",
                "name": "Cluster 9",
                "color": [128, 4, 168]
            },
            {
                "key": "kmeans\tCluster 10",
                "name": "Cluster 10",
                "color": [235, 103, 204]
            }
        ]
    }

    for ci in range(1, 11):
        cluster_cell_barcodes = kmeans_df.index[kmeans_df["Cluster"] == ci].tolist()
        cell_sets["setsTree"][ci]["set"] = cluster_cell_barcodes

    with open(snakemake.output["cells"], "w") as f:
        json.dump(cells, f)
    
    with open(snakemake.output["factors"], "w") as f:
        json.dump(factors, f)
    
    with open(snakemake.output["cell_sets"], "w") as f:
        json.dump(cell_sets, f)
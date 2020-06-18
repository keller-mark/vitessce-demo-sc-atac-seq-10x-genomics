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
    
    # Construct the tree, according to the following schema:
    # https://github.com/hubmapconsortium/vitessce/blob/d5f63aa1d08aa61f6b20f6ad6bbfba5fceb6b5ef/src/schemas/cell_sets.schema.json

    cell_sets = {
        "datatype": "cell",
        "version": "0.1.2",
        "tree": [
            {
                "name": "k-means Clustering",
                "children": [
                    {
                        "name": f"Cluster {cluster_i}",
                        "set": kmeans_df.index[kmeans_df["Cluster"] == cluster_i].tolist()
                    }
                    for cluster_i in range(1, 11)
                ]
            }
        ]
    }

    with open(snakemake.output["cells"], "w") as f:
        json.dump(cells, f)
    
    with open(snakemake.output["cell_sets"], "w") as f:
        json.dump(cell_sets, f)
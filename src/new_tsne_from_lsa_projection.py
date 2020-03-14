import pandas as pd
import numpy as np
from sklearn import manifold


if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0], index_col=0)
    embedding = manifold.TSNE(n_components=2, random_state=0)
    X_transformed = embedding.fit_transform(df.values)

    out_df = pd.DataFrame(index=df.index.values, columns=["TSNE-1", "TSNE-2"], data=X_transformed)
    out_df.index = out_df.index.rename("Barcode")

    out_df.to_csv(snakemake.output[0], index=True)

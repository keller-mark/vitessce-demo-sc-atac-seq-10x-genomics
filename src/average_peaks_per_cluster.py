import collections
import scipy.sparse as sp_sparse
import tables
import pandas as pd
import numpy as np
import h5py

def get_matrix_from_h5(filename):
    # See https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/advanced/h5_matrices
    with h5py.File(filename, "r") as f:
        group = f['matrix']
        feature_group = group['features']
        ids = feature_group['id'][()] # 84626
        #barcodes = group['barcodes'] # 4654
        data = group['data'][()]
        indices = group['indices'][()]
        indptr = group['indptr'][()]
        shape = group['shape'][()] # 84626 x 4654
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        return ids, matrix

def avg_of_cluster(kmeans_df, matrix, ci):
    # The order of rows in the kmeans_df corresponds to the order of rows in the h5 matrix.
    # Obtain the kmeans_df row indices for the cells in the cluster.
    c_cells_i = kmeans_df.index[kmeans_df["Cluster"] == ci].tolist()
    # Select the matrix columns corresponding to the cells in the cluster.
    c_matrix = matrix[:,c_cells_i]
    # Compute the mean for all cells in the matrix
    c_matrix_mean = c_matrix.mean(axis=1) # 84626
    return c_matrix_mean

if __name__ == "__main__":
    kmeans_df = pd.read_csv(snakemake.input["kmeans"])
    ids, matrix = get_matrix_from_h5(snakemake.input["peaks"])
    
    avg_matrix = np.hstack([ avg_of_cluster(kmeans_df, matrix, i) for i in range(1, 11) ])

    # save to output file
    with h5py.File(snakemake.output[0], "w") as f:
        ids_dataset = f.create_dataset("ids", data=ids)
        kmeans_group = f.create_group("kmeans")

        for i in range(avg_matrix.shape[1]):
            cluster_dataset = kmeans_group.create_dataset(str(i+1), data=avg_matrix[:,i])


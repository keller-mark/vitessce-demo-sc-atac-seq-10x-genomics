import collections
import scipy.sparse as sp_sparse
import tables
import pyBigWig
import pandas as pd

FBCMatrix = collections.namedtuple('FBCMatrix', ['ids', 'names', 'barcodes', 'matrix'])
def get_matrix_from_h5(filename):
    # Adapted from https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/advanced/h5_matrices
    with tables.open_file(filename, 'r') as f:
        group = f.get_node(f.root, 'matrix')
        feature_group = getattr(group, 'features')
        ids = getattr(feature_group, 'id') # 84626
        names = getattr(feature_group, 'name') # 84626
        barcodes = getattr(group, 'barcodes') # 4654
        data = getattr(group, 'data')
        indices = getattr(group, 'indices')
        indptr = getattr(group, 'indptr')
        shape = getattr(group, 'shape') # 84626 x 4654
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        # Feature BarCode Matrix format
        return FBCMatrix(ids, names, barcodes, matrix)

def avg_of_cluster(kmeans_df, fbc_matrix, ci):
    # The order of rows in the kmeans_df corresponds to the order of rows in the h5 matrix.
    # Obtain the kmeans_df row indices for the cells in the cluster.
    c_cells_i = kmeans_df.index[kmeans_df["Cluster"] == ci].tolist()
    # Select the matrix columns corresponding to the cells in the cluster.
    c_matrix = fbc_matrix.matrix[:,c_cells_i]
    # Compute the mean for all cells in the matrix
    c_matrix_mean = c_matrix.mean(axis=1) # 84626
    return c_matrix_mean

if __name__ == "__main__":
    kmeans_df = pd.read_csv(snakemake.input["kmeans"])
    fbc_matrix = get_matrix_from_h5(snakemake.input["peaks"])
    avg_matrix = np.hstack([ avg_of_cluster(kmeans_df, fbc_matrix_ci, i) for i in range(1, 11) ])
    
    
    
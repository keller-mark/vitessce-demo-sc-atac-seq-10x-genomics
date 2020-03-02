import h5py
import pyBigWig

if __name__ == "__main__":
    f = h5py.File(snakemake.input[0])
    cluster_i = snakemake.params["i"]

    bw = pyBigWig.open(snakemake.output[0], "w")
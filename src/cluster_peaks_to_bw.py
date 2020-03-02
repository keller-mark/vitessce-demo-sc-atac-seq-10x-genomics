import h5py
import pyBigWig
import numpy as np

if __name__ == "__main__":
    cluster_i = str(snakemake.wildcards["cluster_i"])
    with h5py.File(snakemake.input[0], "r") as f:
        ids = f['ids'][()] # array([b'chr1:237588-237918', b'chr1:565103-565551', ...])
        values = f["kmeans"][cluster_i][:,0]
    
    bw_chroms = [ x[:x.index(b':')].decode("utf-8") for x in ids ]
    bw_starts = [ int(x[x.index(b':')+1:x.index(b'-')]) for x in ids ]
    bw_ends = [ int(x[x.index(b'-')+1:]) for x in ids ]
    bw_values = values

    bw = pyBigWig.open(snakemake.output[0], "w")
    # hg19 header
    bw.addHeader([
        ("chr1", 249250621),
        ("chr2", 243199373),
        ("chr3", 198022430),
        ("chr4", 191154276),
        ("chr5", 180915260),
        ("chr6", 171115067),
        ("chr7", 159138663),
        ("chr8", 146364022),
        ("chr9", 141213431),
        ("chr10", 135534747),
        ("chr11", 135006516),
        ("chr12", 133851895),
        ("chr13", 115169878),
        ("chr14", 107349540),
        ("chr15", 102531392),
        ("chr16", 90354753),
        ("chr17", 81195210),
        ("chr18", 78077248),
        ("chr19", 59128983),
        ("chr20", 63025520),
        ("chr21", 48129895),
        ("chr22", 51304566),
        ("chrX", 155270560),
        ("chrY", 59373566),
        ("chrM", 16571)
    ])
    # values
    bw.addEntries(
        bw_chroms, 
        bw_starts, 
        ends=bw_ends,
        values=bw_values
    )

    bw.close()
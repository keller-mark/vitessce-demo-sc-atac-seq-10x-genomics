import zarr
from numcodecs import Zlib
import bbi
import negspy.coordinates as nc
import numpy as np
import math
import json
from tqdm import tqdm
import resource
import os

def bigwigs_to_zarr(
    input_bigwig_files,
    output_file,
    starting_resolution,
    name
):
    # Short-hand for creating a DirectoryStore with a root group.
    f = zarr.open(output_file, mode='w')
    compressor = Zlib(level=1)

    num_samples = len(input_bigwig_files)

    # Create level zero groups
    chromosomes_group = f.create_group("chromosomes")

    # Prepare to fill in chroms dataset
    chromosomes = nc.get_chromorder('hg38')
    chromosomes = [ str(chr_name) for chr_name in chromosomes[:25] ] # TODO: should more than chr1-chrM be used?
    num_chromosomes = len(chromosomes)
    chroms_length_arr = np.array([ nc.get_chrominfo('hg38').chrom_lengths[x] for x in chromosomes ], dtype="i8")
    chroms_cumsum_arr = np.concatenate((np.array([0]), np.cumsum(chroms_length_arr)))

    chromosomes_set = set(chromosomes)
    chrom_name_to_length = dict(zip(chromosomes, chroms_length_arr))
    chrom_name_to_cumsum = dict(zip(chromosomes, chroms_cumsum_arr))

    
    # Prepare to fill in resolutions dataset
    resolutions = [ starting_resolution*(2**x) for x in range(16)]
        
    # Create each chromosome dataset.
    for chr_name, chr_len in chrom_name_to_length.items():
        chr_group = chromosomes_group.create_group(chr_name)
        # Create each resolution group.
        for resolution in resolutions:
            chr_shape = (num_samples, math.ceil(chr_len / resolution))
            chr_group.create_dataset(str(resolution), shape=chr_shape, dtype="f4", fill_value=np.nan, compressor=compressor)
    
    # Fill in data for each bigwig file.
    for bw_index, bw_file in tqdm(list(enumerate(input_bigwig_files)), desc='bigwigs'):
        if bbi.is_bigwig(bw_file):
            chromsizes = bbi.chromsizes(bw_file)
            matching_chromosomes = set(chromsizes.keys()).intersection(chromosomes_set)

            # Fill in data for each resolution of a bigwig file.
            for resolution in resolutions:
                # Fill in data for each chromosome of a resolution of a bigwig file.
                for chr_name in matching_chromosomes:
                    chr_len = chrom_name_to_length[chr_name]
                    chr_shape = (num_samples, math.ceil(chr_len / resolution))
                    arr = bbi.fetch(bw_file, chr_name, 0, chr_len, chr_shape[1], summary="sum")
                    chromosomes_group[chr_name][str(resolution)][bw_index,:] = arr
        else:
            print(f"{bw_file} not is_bigwig")
    
    max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(max_mem)
    
    # Append metadata to the top resolution row_infos attribute.
    row_infos = []
    for bw_index, bw_file in enumerate(input_bigwig_files):
        row_infos.append({ "cluster": int(bw_index+1), "file": os.path.basename(bw_file) })
    
    # f.attrs should contain all tileset_info properties
    # For zarr, more attributes are used here to allow "serverless"
    f.attrs['row_infos'] = row_infos
    f.attrs['resolutions'] = sorted(resolutions, reverse=True)
    f.attrs['shape'] = [ num_samples, 256 ]
    f.attrs['name'] = name
    f.attrs['coordSystem'] = "hg38"
    
    # https://github.com/zarr-developers/zarr-specs/issues/50
    f.attrs['multiscales'] = [
        {
            "version": "0.1",
            "name": chr_name,
            "datasets": [
                { "path": f"chromosomes/{chr_name}/{resolution}" }
                for resolution in sorted(resolutions, reverse=True)
            ],
            "type": "zarr-multivec",
            "metadata": {
                "chromoffset": int(chrom_name_to_cumsum[chr_name]),
                "chromsize": int(chr_len),
            }
        }
        for (chr_name, chr_len) in list(zip(chromosomes, chroms_length_arr))
    ]


if __name__ == "__main__":
    bigwigs_to_zarr(
        list(snakemake.input),
        snakemake.output[0],
        25,
        "k-means clusters"
    )
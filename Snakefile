from os.path import join

DATA_DIR = "data"
RAW_DIR = join(DATA_DIR, "raw")
PROCESSED_DIR = join(DATA_DIR, "processed")
SRC_DIR = "src"

BASE_URL = "http://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_pbmc_5k"

rule all:
    input:
        join(PROCESSED_DIR, "pbmc_10x.cells.json"),
        join(PROCESSED_DIR, "pbmc_10x.factors.json"),
        join(PROCESSED_DIR, "pbmc_10x.cell_sets.json"),
        [ join(PROCESSED_DIR, f"pbmc_10x_peaks_{i}.bw") for i in range(1, 11) ]

rule cluster_peaks_to_bw:
    input:
        join(RAW_DIR, "pbmc_10x_peaks_avg_per_cluster.h5")
    output:
        join(PROCESSED_DIR, "pbmc_10x_peaks_{cluster_i}.bw")
    script:
        join(SRC_DIR, "cluster_peaks_to_bw.py")

rule average_peaks_per_cluster:
    input:
        peaks=join(RAW_DIR, "filtered_peak_bc_matrix.h5"),
        kmeans=join(RAW_DIR, "analysis", "clustering", "kmeans_10_clusters", "clusters.csv")
    output:
        join(RAW_DIR, "pbmc_10x_peaks_avg_per_cluster.h5")
    script:
        join(SRC_DIR, "average_peaks_per_cluster.py")

rule process_cells_and_factors:
    input:
        tsne_1=join(RAW_DIR, "analysis", "tsne", "2_components", "projection.csv"),
        tsne_2=join(RAW_DIR, "tsne_from_lsa.csv"),
        kmeans=join(RAW_DIR, "analysis", "clustering", "kmeans_10_clusters", "clusters.csv")
    output:
        cells=join(PROCESSED_DIR, "pbmc_10x.cells.json"),
        factors=join(PROCESSED_DIR, "pbmc_10x.factors.json"),
        cell_sets=join(PROCESSED_DIR, "pbmc_10x.cell_sets.json")
    script:
        join(SRC_DIR, "process_cells_and_factors.py")

rule new_tsne_from_lsa_projection:
    input:
        join(RAW_DIR, "analysis", "lsa", "15_components", "projection.csv")
    output:
        join(RAW_DIR, "tsne_from_lsa.csv")
    script:
        join(SRC_DIR, "new_tsne_from_lsa_projection.py")

rule untar:
    input:
        analysis_tar=join(RAW_DIR, "pbmc_analysis.tar.gz")
    output:
        join(RAW_DIR, "analysis", "tsne", "2_components", "projection.csv"),
        join(RAW_DIR, "analysis", "clustering", "kmeans_10_clusters", "clusters.csv"),
        join(RAW_DIR, "analysis", "lsa", "15_components", "projection.csv")
    shell:
        """
        tar -xvzf {input.analysis_tar} -C {RAW_DIR}
        """

rule download:
    output:
        analysis_tar=join(RAW_DIR, "pbmc_analysis.tar.gz"),
        f_peak_bc_matrix_h5=join(RAW_DIR, "filtered_peak_bc_matrix.h5"),
        peak_annotation_tsv=join(RAW_DIR, "peak_annotation.tsv"),
        peaks_bed=join(RAW_DIR, "peaks.bed"),
        singlecell_csv=join(RAW_DIR, "singlecell.csv"),
        summary_csv=join(RAW_DIR, "summary.csv")
    shell:
        """
        curl -L -o {output.analysis_tar} {BASE_URL}/atac_v1_pbmc_5k_analysis.tar.gz && \
        curl -L -o {output.f_peak_bc_matrix_h5} {BASE_URL}/atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5 && \
        curl -L -o {output.peak_annotation_tsv} {BASE_URL}/atac_v1_pbmc_5k_peak_annotation.tsv && \
        curl -L -o {output.peaks_bed} {BASE_URL}/atac_v1_pbmc_5k_peaks.bed && \
        curl -L -o {output.singlecell_csv} {BASE_URL}/atac_v1_pbmc_5k_singlecell.csv && \
        curl -L -o {output.summary_csv} {BASE_URL}/atac_v1_pbmc_5k_summary.csv
        """
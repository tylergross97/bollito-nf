process SEURAT_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::r-seurat=5.1.0 conda-forge::r-dplyr conda-forge::r-data.table conda-forge::r-ggplot2'

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("${meta.id}_seurat_post-qc-filtered.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("dplyr"))

    message("=== SEURAT GENE FILTER: ${meta.id} ===")

    random_seed <- ${params.random_seed}
    set.seed(random_seed)

    gene_filter    <- strsplit("${params.filter_genes}", ",")[[1]]
    filter_out     <- ${params.filter_out ? 'TRUE' : 'FALSE'}
    threshold      <- ${params.filter_threshold}

    seurat_obj <- readRDS("${rds}")

    # For each gene, check if cells express it above threshold
    genes_present <- intersect(gene_filter, rownames(seurat_obj))
    if (length(genes_present) == 0) {
        message("None of the filter genes found in the dataset. Skipping filter.")
    } else {
        # Get expression for filter genes
        expr_matrix <- GetAssayData(seurat_obj, layer = "counts")[genes_present, , drop = FALSE]
        cells_expressing <- colSums(expr_matrix > 0) > 0
        cells_above_thresh <- colSums(expr_matrix) >= threshold

        if (filter_out) {
            # Keep cells WITH expression >= threshold
            seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_above_thresh])
        } else {
            # Remove cells WITH expression >= threshold (keep those below)
            seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[!cells_above_thresh])
        }
    }

    saveRDS(seurat_obj, file = "${meta.id}_seurat_post-qc-filtered.rds")
    message("Gene filtering complete for ${meta.id}")
    """
}

process SEURAT_DEGS {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::r-seurat=5.1.0 conda-forge::r-dplyr conda-forge::r-data.table conda-forge::r-ggplot2 bioconda::bioconductor-biocparallel conda-forge::r-openxlsx conda-forge::r-future'

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("${meta.id}_seurat_degs.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("BiocParallel"))
    suppressMessages(library("openxlsx"))
    suppressMessages(library("future"))

    message("=== SEURAT DEGs: ${meta.id} ===")

    random_seed   <- ${params.random_seed}
    set.seed(random_seed)

    selected_cond <- strsplit("${params.selected_cond}", ",")[[1]]
    deg_ranking   <- ${params.deg_ranking ? 'TRUE' : 'FALSE'}
    deg_test      <- "${params.deg_test}"

    plan("multicore", workers = ${task.cpus})
    options(future.globals.maxSize = ${task.memory.toMega()} * 1024^2)

    seurat_obj <- readRDS("${rds}")

    # Parse selected condition/resolution
    # selected_cond can be a resolution (e.g., "0.2") or a metadata column (e.g., "condition")
    for (cond in selected_cond) {
        # Try to use it as a resolution
        res_col <- paste0("RNA_snn_res.", cond)
        if (!res_col %in% colnames(seurat_obj@meta.data)) {
            res_col <- paste0("SCT_snn_res.", cond)
        }
        if (res_col %in% colnames(seurat_obj@meta.data)) {
            Idents(seurat_obj) <- res_col
        } else if (cond %in% colnames(seurat_obj@meta.data)) {
            Idents(seurat_obj) <- cond
        } else {
            message("Condition/resolution not found: ", cond, ". Skipping.")
            next
        }

        # Find all markers
        tryCatch({
            markers <- FindAllMarkers(seurat_obj, test.use = deg_test, verbose = FALSE)
            if (nrow(markers) > 0) {
                write_xlsx(markers, path = paste0("${meta.id}_markers_", cond, ".xlsx"))
            }
        }, error = function(e) {
            message("FindAllMarkers failed for ", cond, ": ", e\$message)
        })

        # Gene ranking (for GSEA)
        if (deg_ranking) {
            tryCatch({
                markers_all <- FindAllMarkers(seurat_obj, test.use = deg_test,
                                              min.pct = 0, logfc.threshold = 0, verbose = FALSE)
                if (nrow(markers_all) > 0) {
                    write_xlsx(markers_all, path = paste0("${meta.id}_ranking_", cond, ".xlsx"))
                }
            }, error = function(e) {
                message("Gene ranking failed for ", cond, ": ", e\$message)
            })
        }
    }

    saveRDS(seurat_obj, file = "${meta.id}_seurat_degs.rds")
    message("DEG analysis complete for ${meta.id}")
    """
}

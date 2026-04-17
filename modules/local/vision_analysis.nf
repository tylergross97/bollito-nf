process VISION_ANALYSIS {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::r-seurat=5.1.0 conda-forge::r-ggplot2 conda-forge::r-viridis conda-forge::r-rcolorbrewer conda-forge::r-scales conda-forge::r-devtools'

    input:
    tuple val(meta), path(rds)
    path gmt

    output:
    tuple val(meta), path("${meta.id}_vision_object.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("viridis"))
    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("scales"))

    # VISION might need to be installed from GitHub
    tryCatch({
        suppressMessages(library("VISION"))
    }, error = function(e) {
        message("VISION library not available. Attempting install...")
        devtools::install_github("YosefLab/VISION", upgrade = "never")
        suppressMessages(library("VISION"))
    })

    message("=== VISION FUNCTIONAL ANALYSIS: ${meta.id} ===")

    random_seed     <- ${params.random_seed}
    set.seed(random_seed)

    meta_columns    <- strsplit("${params.vision_meta_columns}", ",")[[1]]
    selected_res    <- ${params.vision_res}
    regress_out     <- ${params.regress_out ? 'TRUE' : 'FALSE'}
    vars_to_regress <- strsplit("${params.vars_to_regress}", ",")[[1]]
    regress_cc      <- ${params.regress_cell_cycle ? 'TRUE' : 'FALSE'}

    options(future.globals.maxSize = ${task.memory.toMega()} * 1024^2)

    seurat_obj <- readRDS("${rds}")

    # Set resolution
    res_col <- paste0("RNA_snn_res.", selected_res)
    if (!res_col %in% colnames(seurat_obj@meta.data)) {
        res_col <- paste0("SCT_snn_res.", selected_res)
    }

    # Build metadata for VISION
    meta_df <- seurat_obj@meta.data
    cols_to_keep <- intersect(c(meta_columns, res_col), colnames(meta_df))
    meta_for_vision <- meta_df[, cols_to_keep, drop = FALSE]

    # Get expression matrix
    expr <- GetAssayData(seurat_obj, layer = "data")

    # Create VISION object
    vis <- Vision(expr,
                  signatures = "${gmt}",
                  meta = meta_for_vision)

    # Run analysis
    vis <- analyze(vis)

    saveRDS(vis, file = "${meta.id}_vision_object.rds")
    message("VISION analysis complete for ${meta.id}")
    """
}

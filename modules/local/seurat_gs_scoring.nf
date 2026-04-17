process SEURAT_GS_SCORING {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::r-seurat=5.1.0 conda-forge::r-dplyr conda-forge::r-data.table conda-forge::r-ggplot2 conda-forge::r-clustree conda-forge::r-patchwork bioconda::r-qusage=2.38.0'

    input:
    tuple val(meta), path(rds)
    path gmt

    output:
    tuple val(meta), path("${meta.id}_seurat_complete.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("qusage"))
    suppressMessages(library("patchwork"))

    message("=== SEURAT GENE SET SCORING: ${meta.id} ===")

    random_seed <- ${params.random_seed}
    set.seed(random_seed)

    resolutions        <- as.numeric(strsplit("${params.resolutions}", ",")[[1]])
    norm_type          <- "${params.norm_type}"
    geneset_percentage <- ${params.geneset_percentage}

    options(future.globals.maxSize = ${task.memory.toMega()} * 1024^2)

    seurat_obj <- readRDS("${rds}")

    # Read GMT
    genesets <- read.gmt("${gmt}")

    # Filter genesets by expression coverage
    all_genes <- rownames(seurat_obj)
    genesets_filtered <- list()
    for (gs_name in names(genesets)) {
        gs_genes <- genesets[[gs_name]]
        pct_expressed <- sum(gs_genes %in% all_genes) / length(gs_genes)
        if (pct_expressed >= geneset_percentage) {
            genesets_filtered[[gs_name]] <- gs_genes
        }
    }

    message(paste0("Gene sets passing filter: ", length(genesets_filtered), "/", length(genesets)))

    # Score each geneset
    if (length(genesets_filtered) > 0) {
        for (gs_name in names(genesets_filtered)) {
            gs_genes <- intersect(genesets_filtered[[gs_name]], all_genes)
            if (length(gs_genes) >= 2) {
                tryCatch({
                    seurat_obj <- AddModuleScore(seurat_obj,
                                                 features = list(gs_genes),
                                                 name = gs_name,
                                                 seed = random_seed)
                }, error = function(e) {
                    message("Scoring failed for ", gs_name, ": ", e\$message)
                })
            }
        }
    }

    saveRDS(seurat_obj, file = "${meta.id}_seurat_complete.rds")
    message("Gene set scoring complete for ${meta.id}")
    """
}

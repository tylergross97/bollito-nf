process SEURAT_POSTQC {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::r-seurat=5.1.0 conda-forge::r-dplyr conda-forge::r-data.table conda-forge::r-ggplot2 conda-forge::r-patchwork'

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("${meta.id}_seurat_post-qc.rds"), emit: rds

    script:
    def min_count_r = params.min_count ? "${params.min_count}" : "NULL"
    def max_count_r = params.max_count ? "${params.max_count}" : "NULL"
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("patchwork"))

    message("=== SEURAT POST-QC: ${meta.id} ===")

    random_seed <- ${params.random_seed}
    set.seed(random_seed)

    min_feat   <- ${params.min_feat}
    max_feat   <- ${params.max_feat}
    min_count  <- ${min_count_r}
    max_count  <- ${max_count_r}
    mit_pct    <- ${params.mit_pct}
    ribo_pct   <- ${params.ribo_pct}
    write_table <- ${params.write_table ? 'TRUE' : 'FALSE'}

    # Load pre-QC object
    seurat_obj <- readRDS("${rds}")

    # Count cells before filtering
    ncells_pre <- ncol(seurat_obj)

    # Apply filters
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= min_feat & nFeature_RNA <= max_feat)
    seurat_obj <- subset(seurat_obj, subset = percent.mt <= mit_pct)
    seurat_obj <- subset(seurat_obj, subset = percent.ribo <= ribo_pct)

    if (!is.null(min_count)) {
        seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= min_count)
    }
    if (!is.null(max_count)) {
        seurat_obj <- subset(seurat_obj, subset = nCount_RNA <= max_count)
    }

    ncells_post <- ncol(seurat_obj)
    message(paste0("Cells: ", ncells_pre, " -> ", ncells_post, " (removed ", ncells_pre - ncells_post, ")"))

    # Post-filter QC plots (direct ggplot2 to avoid Seurat v5 VlnPlot/FetchData API issues)
    qc_meta <- seurat_obj@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")]
    vln_plots <- lapply(colnames(qc_meta), function(feat) {
        df <- data.frame(x = "cells", y = qc_meta[[feat]])
        ggplot(df, aes(x = x, y = y)) +
            geom_violin(fill = "grey80") +
            geom_jitter(size = 0.1, alpha = 0.3, width = 0.2) +
            labs(title = feat, x = "", y = "") +
            theme_classic() +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    })
    pdf("${meta.id}_vlnplot_QC_postfilt.pdf", width = 12, height = 6)
    print(patchwork::wrap_plots(vln_plots, ncol = 4))
    dev.off()

    # Stats table
    if (write_table) {
        stats <- data.frame(
            metric = c("cells_pre_filter", "cells_post_filter", "cells_removed",
                        "min_feat", "max_feat", "mit_pct", "ribo_pct"),
            value  = c(ncells_pre, ncells_post, ncells_pre - ncells_post,
                        min_feat, max_feat, mit_pct, ribo_pct)
        )
        write.table(stats, file = "${meta.id}_pre_vs_post_stats.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)
    }

    saveRDS(seurat_obj, file = "${meta.id}_seurat_post-qc.rds")
    message("Post-QC filtering complete for ${meta.id}")
    """
}

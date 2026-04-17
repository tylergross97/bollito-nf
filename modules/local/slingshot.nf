process SLINGSHOT {
    tag "${meta.id}"
    label 'process_medium'

    conda 'bioconda::bioconductor-slingshot=2.12.0 bioconda::bioconductor-singlecellexperiment conda-forge::r-seurat=5.1.0 conda-forge::r-rcolorbrewer conda-forge::r-gam conda-forge::r-pheatmap bioconda::bioconductor-genomeinfodb'

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("${meta.id}_slingshot_sce_objects.RData"), emit: rdata

    script:
    def start_clus_r = params.start_clus ? "\"${params.start_clus}\"" : "NULL"
    def end_clus_r   = params.end_clus   ? "\"${params.end_clus}\""   : "NULL"
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("slingshot"))
    suppressMessages(library("SingleCellExperiment"))
    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("gam"))
    suppressMessages(library("pheatmap"))

    message("=== SLINGSHOT TRAJECTORY: ${meta.id} ===")

    random_seed    <- ${params.random_seed}
    set.seed(random_seed)

    selected_res   <- ${params.slingshot_res}
    start_clus     <- ${start_clus_r}
    end_clus       <- ${end_clus_r}
    n_var_genes    <- ${params.n_var_genes}
    n_plotted_genes <- ${params.n_plotted_genes}
    pc             <- ${params.principal_components}
    graphics       <- ${params.graphics ? 'TRUE' : 'FALSE'}

    options(future.globals.maxSize = ${task.memory.toMega()} * 1024^2)

    seurat_obj <- readRDS("${rds}")

    # Set identity to selected resolution
    res_col <- paste0("RNA_snn_res.", selected_res)
    if (!res_col %in% colnames(seurat_obj@meta.data)) {
        res_col <- paste0("SCT_snn_res.", selected_res)
    }
    Idents(seurat_obj) <- res_col

    # Convert to SCE
    sce <- as.SingleCellExperiment(seurat_obj)

    # Run slingshot
    slingshot_args <- list(data = sce, clusterLabels = "ident", reducedDim = "PCA")
    if (!is.null(start_clus)) slingshot_args\$start.clus <- start_clus
    if (!is.null(end_clus))   slingshot_args\$end.clus   <- end_clus

    sce <- do.call(slingshot, slingshot_args)

    # Plot trajectories
    colors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(seurat_obj@active.ident)))

    pdf("${meta.id}_slingshot_trajectory.pdf", width = 10, height = 8)
    plot(reducedDims(sce)\$UMAP, col = colors[sce\$ident], pch = 16, cex = 0.5)
    lines(SlingshotDataSet(sce), lwd = 2, col = "black")
    dev.off()

    # Save
    save(sce, seurat_obj, file = "${meta.id}_slingshot_sce_objects.RData")
    message("Slingshot trajectory inference complete for ${meta.id}")
    """
}

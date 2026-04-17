process RNA_VELOCITY {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::r-seurat=5.1.0 bioconda::r-velocyto.r=0.6 conda-forge::r-rcolorbrewer conda-forge::r-future bioconda::r-seuratdisk=0.0.0.9021'

    input:
    tuple val(meta), path(rds)
    tuple val(meta2), path(velocyto_dir)

    output:
    tuple val(meta), path("${meta.id}_seurat_velocity.rds"), emit: rds

    script:
    def n_cells_r = params.downsample_n_cells ? "${params.downsample_n_cells}" : "NULL"
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("velocyto.R"))
    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("future"))

    # Try to load SeuratWrappers
    tryCatch({
        suppressMessages(library("SeuratWrappers"))
    }, error = function(e) {
        message("SeuratWrappers not available, installing...")
        remotes::install_github("satijalab/seurat-wrappers", upgrade = "never")
        suppressMessages(library("SeuratWrappers"))
    })

    message("=== RNA VELOCITY: ${meta.id} ===")

    random_seed  <- ${params.random_seed}
    set.seed(random_seed)

    selected_res <- ${params.velocyto_res}
    downsampling <- ${params.downsampling ? 'TRUE' : 'FALSE'}
    n_cells      <- ${n_cells_r}

    plan("multicore", workers = ${task.cpus})
    options(future.globals.maxSize = ${task.memory.toMega()} * 1024^2)

    seurat_obj <- readRDS("${rds}")

    # Set resolution
    res_col <- paste0("RNA_snn_res.", selected_res)
    if (!res_col %in% colnames(seurat_obj@meta.data)) {
        res_col <- paste0("SCT_snn_res.", selected_res)
    }
    Idents(seurat_obj) <- res_col

    # Load velocyto matrices
    spliced_dir   <- file.path("${velocyto_dir}", "spliced")
    unspliced_dir <- file.path("${velocyto_dir}", "unspliced")
    ambiguous_dir <- file.path("${velocyto_dir}", "ambiguous")

    spliced   <- ReadMtx(mtx = file.path(spliced_dir, "matrix.mtx"),
                         cells = file.path(spliced_dir, "barcodes.tsv"),
                         features = file.path(spliced_dir, "genes.tsv"))
    unspliced <- ReadMtx(mtx = file.path(unspliced_dir, "matrix.mtx"),
                         cells = file.path(unspliced_dir, "barcodes.tsv"),
                         features = file.path(unspliced_dir, "genes.tsv"))

    # Downsample if requested
    if (downsampling && !is.null(n_cells) && n_cells < ncol(seurat_obj)) {
        cells_to_keep <- sample(colnames(seurat_obj), n_cells)
        seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    }

    # Add spliced/unspliced assays
    common_cells <- intersect(colnames(seurat_obj), colnames(spliced))
    common_genes <- intersect(rownames(seurat_obj), rownames(spliced))

    seurat_obj[["spliced"]]   <- CreateAssayObject(counts = spliced[common_genes, common_cells])
    seurat_obj[["unspliced"]] <- CreateAssayObject(counts = unspliced[common_genes, common_cells])

    # Run velocity estimation
    tryCatch({
        seurat_obj <- RunVelocity(seurat_obj,
                                   deltaT = 1,
                                   kCells = 25,
                                   fit.quantile = 0.02,
                                   spliced.average = "spliced",
                                   unspliced.average = "unspliced")

        colors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(levels(seurat_obj)))

        pdf("${meta.id}_velocity_umap.pdf", width = 10, height = 8)
        show.velocity.on.embedding.cor(
            emb = Embeddings(seurat_obj, reduction = "umap"),
            vel = Tool(seurat_obj, slot = "RunVelocity"),
            n = 200, scale = "sqrt",
            cell.colors = ac(x = colors[as.numeric(Idents(seurat_obj))], alpha = 0.5),
            cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE,
            min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1
        )
        dev.off()
    }, error = function(e) {
        message("Velocity estimation failed: ", e\$message)
    })

    saveRDS(seurat_obj, file = "${meta.id}_seurat_velocity.rds")
    message("RNA velocity complete for ${meta.id}")
    """
}

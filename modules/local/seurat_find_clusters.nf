process SEURAT_FIND_CLUSTERS {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::r-seurat=5.1.0 conda-forge::r-dplyr conda-forge::r-data.table conda-forge::r-ggplot2 conda-forge::r-clustree conda-forge::r-cluster conda-forge::r-writexl conda-forge::r-future conda-forge::r-patchwork conda-forge::r-lisi bioconda::r-seuratdisk=0.0.0.9021'

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("${meta.id}_seurat_find-clusters.rds"),  emit: rds
    tuple val(meta), path("${meta.id}_seurat_find-clusters.h5ad"), emit: h5ad

    script:
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("SeuratDisk"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("clustree"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("cluster"))
    suppressMessages(library("writexl"))
    suppressMessages(library("future"))
    suppressMessages(library("patchwork"))
    suppressMessages(library("lisi"))

    message("=== SEURAT FIND CLUSTERS: ${meta.id} ===")

    random_seed <- ${params.random_seed}
    set.seed(random_seed)

    pc            <- ${params.principal_components}
    resolutions   <- as.numeric(strsplit("${params.resolutions}", ",")[[1]])
    k_neighbors   <- ${params.k_neighbors}
    batch_metadata <- "${params.batch_metadata}"

    plan("multicore", workers = ${task.cpus})
    options(future.globals.maxSize = ${task.memory.toMega()} * 1024^2)

    seurat_obj <- readRDS("${rds}")

    # Find neighbors
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pc, k.param = k_neighbors, verbose = FALSE)

    # Find clusters at each resolution
    for (res in resolutions) {
        seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    }

    # Run UMAP
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:pc, verbose = FALSE)

    # Clustree plot
    pdf("${meta.id}_clustree.pdf", width = 12, height = 10)
    tryCatch({
        print(clustree(seurat_obj, prefix = "RNA_snn_res."))
    }, error = function(e) {
        tryCatch({
            print(clustree(seurat_obj, prefix = "SCT_snn_res."))
        }, error = function(e2) {
            message("Clustree failed: ", e2\$message)
        })
    })
    dev.off()

    # UMAP plots for each resolution
    for (res in resolutions) {
        res_col <- paste0("RNA_snn_res.", res)
        if (!res_col %in% colnames(seurat_obj@meta.data)) {
            res_col <- paste0("SCT_snn_res.", res)
        }
        if (res_col %in% colnames(seurat_obj@meta.data)) {
            pdf(paste0("${meta.id}_umap_res_", res, ".pdf"), width = 8, height = 6)
            Idents(seurat_obj) <- res_col
            print(DimPlot(seurat_obj, reduction = "umap", label = TRUE))
            dev.off()
        }
    }

    # Compute LISI if batch metadata exists
    tryCatch({
        if (batch_metadata %in% colnames(seurat_obj@meta.data)) {
            lisi_res <- compute_lisi(Embeddings(seurat_obj, "umap"),
                                     seurat_obj@meta.data,
                                     batch_metadata)
            seurat_obj\$LISI <- lisi_res[, 1]
        }
    }, error = function(e) {
        message("LISI computation failed: ", e\$message)
    })

    # Save RDS
    saveRDS(seurat_obj, file = "${meta.id}_seurat_find-clusters.rds")

    # Save h5ad for scanpy interoperability
    tryCatch({
        SaveH5Seurat(seurat_obj, filename = "${meta.id}_seurat_find-clusters.h5seurat", overwrite = TRUE)
        Convert("${meta.id}_seurat_find-clusters.h5seurat", dest = "h5ad", overwrite = TRUE)
    }, error = function(e) {
        message("h5ad conversion failed: ", e\$message)
        # Create empty placeholder
        file.create("${meta.id}_seurat_find-clusters.h5ad")
    })

    message("Clustering complete for ${meta.id}")
    """
}

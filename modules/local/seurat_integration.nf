process SEURAT_INTEGRATION {
    tag 'integrated'
    label 'process_high'

    conda 'conda-forge::r-seurat=5.1.0 conda-forge::r-ggplot2 conda-forge::r-stringr conda-forge::r-future conda-forge::r-sctransform'

    input:
    path rds_files     // collected list of per-sample RDS files
    val sample_ids     // collected list of sample IDs
    val meta           // meta map passed from workflow

    output:
    tuple val(meta), path("integrated_seurat_normalized-pcs.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("stringr"))
    suppressMessages(library("future"))

    message("=== SEURAT INTEGRATION ===")

    random_seed     <- ${params.random_seed}
    set.seed(random_seed)
    norm_type       <- "${params.norm_type}"
    vf              <- ${params.variable_features ? 'TRUE' : 'FALSE'}
    vars_to_regress <- strsplit("${params.vars_to_regress}", ",")[[1]]
    gene_case       <- "${params.gene_case}"
    write_table     <- ${params.write_table ? 'TRUE' : 'FALSE'}

    plan("multicore", workers = ${task.cpus})
    options(future.globals.maxSize = ${task.memory.toMega()} * 1024^2)

    # Load all sample objects
    rds_paths    <- strsplit("${rds_files.join(',')}", ",")[[1]]
    sample_names <- strsplit("${sample_ids.join(',')}", ",")[[1]]

    objects <- list()
    for (i in seq_along(rds_paths)) {
        obj <- readRDS(trimws(rds_paths[i]))
        obj\$orig.sample <- sample_names[i]
        objects[[sample_names[i]]] <- obj
    }

    if (length(objects) < 2) {
        stop("Integration requires at least 2 samples. Found: ", length(objects))
    }

    # Normalize each object individually
    if (norm_type == "SCT") {
        objects <- lapply(objects, function(obj) {
            SCTransform(obj, vars.to.regress = vars_to_regress, verbose = FALSE)
        })
        features <- SelectIntegrationFeatures(object.list = objects, nfeatures = 3000)
        objects  <- PrepSCTIntegration(object.list = objects, anchor.features = features)
        anchors  <- FindIntegrationAnchors(object.list = objects, normalization.method = "SCT",
                                           anchor.features = features)
        integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
    } else {
        objects <- lapply(objects, function(obj) {
            obj <- NormalizeData(obj, verbose = FALSE)
            obj <- FindVariableFeatures(obj, verbose = FALSE)
            return(obj)
        })
        features <- SelectIntegrationFeatures(object.list = objects)
        anchors  <- FindIntegrationAnchors(object.list = objects, anchor.features = features)
        integrated <- IntegrateData(anchorset = anchors)
        integrated <- ScaleData(integrated, verbose = FALSE)
    }

    # PCA
    if (vf) {
        integrated <- RunPCA(integrated, verbose = FALSE)
    } else {
        integrated <- RunPCA(integrated, features = rownames(integrated), verbose = FALSE)
    }

    # Plots
    pdf("integrated_elbowplot.pdf", width = 8, height = 5)
    print(ElbowPlot(integrated))
    dev.off()

    pdf("integrated_dimplot.pdf", width = 8, height = 6)
    print(DimPlot(integrated, reduction = "pca", group.by = "orig.sample"))
    dev.off()

    saveRDS(integrated, file = "integrated_seurat_normalized-pcs.rds")
    message("Integration complete")
    """
}

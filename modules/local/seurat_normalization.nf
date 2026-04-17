process SEURAT_NORMALIZATION {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::r-seurat=5.1.0 conda-forge::r-dplyr conda-forge::r-data.table conda-forge::r-ggplot2 conda-forge::r-stringr conda-forge::r-future conda-forge::r-sctransform'

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("${meta.id}_seurat_normalized-pcs.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("stringr"))
    suppressMessages(library("future"))

    message("=== SEURAT NORMALIZATION: ${meta.id} ===")

    random_seed    <- ${params.random_seed}
    set.seed(random_seed)

    norm_type      <- "${params.norm_type}"
    vf             <- ${params.variable_features ? 'TRUE' : 'FALSE'}
    regress_out    <- ${params.regress_out ? 'TRUE' : 'FALSE'}
    vars_to_regress <- strsplit("${params.vars_to_regress}", ",")[[1]]
    regress_cc     <- ${params.regress_cell_cycle ? 'TRUE' : 'FALSE'}
    regress_merge  <- ${params.regress_merge_effect ? 'TRUE' : 'FALSE'}
    gene_case      <- "${params.gene_case}"
    write_table    <- ${params.write_table ? 'TRUE' : 'FALSE'}

    plan("multicore", workers = ${task.cpus})
    options(future.globals.maxSize = ${task.memory.toMega()} * 1024^2)

    seurat_obj <- readRDS("${rds}")

    # Normalization
    if (norm_type == "SCT") {
        if (regress_out) {
            seurat_obj <- SCTransform(seurat_obj, vars.to.regress = vars_to_regress, verbose = FALSE)
        } else {
            seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
        }
    } else {
        seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
        seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
        if (regress_out) {
            seurat_obj <- ScaleData(seurat_obj, vars.to.regress = vars_to_regress, verbose = FALSE)
        } else {
            seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
        }
    }

    # Cell cycle scoring
    if (gene_case == "lowercase") {
        s.genes   <- tolower(cc.genes\$s.genes)
        g2m.genes <- tolower(cc.genes\$g2m.genes)
    } else if (gene_case == "titlecase") {
        s.genes   <- str_to_title(cc.genes\$s.genes)
        g2m.genes <- str_to_title(cc.genes\$g2m.genes)
    } else {
        s.genes   <- cc.genes\$s.genes
        g2m.genes <- cc.genes\$g2m.genes
    }
    tryCatch({
        seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
    }, error = function(e) {
        message("Cell cycle scoring failed: ", e\$message)
    })

    # PCA
    if (vf) {
        seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    } else {
        seurat_obj <- RunPCA(seurat_obj, features = rownames(seurat_obj), verbose = FALSE)
    }

    # Plots
    pdf("${meta.id}_elbowplot.pdf", width = 8, height = 5)
    print(ElbowPlot(seurat_obj))
    dev.off()

    pdf("${meta.id}_dimplot.pdf", width = 8, height = 6)
    print(DimPlot(seurat_obj, reduction = "pca"))
    dev.off()

    pdf("${meta.id}_cell_cycle_dimplot.pdf", width = 8, height = 6)
    tryCatch({
        print(DimPlot(seurat_obj, reduction = "pca", group.by = "Phase"))
    }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Cell cycle data not available")
    })
    dev.off()

    saveRDS(seurat_obj, file = "${meta.id}_seurat_normalized-pcs.rds")
    message("Normalization complete for ${meta.id}")
    """
}

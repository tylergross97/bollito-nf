process SEURAT_QC {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::r-seurat=5.1.0 conda-forge::r-dplyr conda-forge::r-data.table conda-forge::r-ggplot2 conda-forge::r-stringr conda-forge::r-matrix conda-forge::r-patchwork conda-forge::r-reticulate'

    input:
    tuple val(meta), path(input_data)
    path samples_file

    output:
    tuple val(meta), path("${meta.id}_seurat_pre-qc.rds"), emit: rds

    script:
    def input_dir = params.input_type == 'fastq' ? "${input_data}/Gene/raw" : "${input_data}"
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("data.table"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("stringr"))
    suppressMessages(library("Matrix"))
    suppressMessages(library("patchwork"))

    message("=== SEURAT QC: ${meta.id} ===")

    # Parameters from Nextflow
    data_dir       <- "${input_dir}"
    sample_id      <- "${meta.id}"
    input_type     <- "${params.input_type}"
    technology     <- "${params.technology}"
    random_seed    <- ${params.random_seed}
    gene_case      <- "${params.gene_case}"
    min_cells      <- ${params.min_cells_per_gene}

    set.seed(random_seed)

    # Set regexp for QC variables
    if (gene_case == "lowercase") {
        mito_grep <- "^mt-"
        ribo_grep <- "^rp[sl][[:digit:]]"
    } else if (gene_case == "titlecase") {
        mito_grep <- "^mt-"
        ribo_grep <- "^Rp[sl][[:digit:]]"
    } else {
        mito_grep <- "^MT-"
        ribo_grep <- "^RP[SL][[:digit:]]"
    }

    # Load data
    if (input_type == "fastq" || technology == "10x") {
        expression_matrix <- Read10X(data.dir = data_dir)
        if (is.list(expression_matrix)) {
            expression_matrix <- expression_matrix[["Gene Expression"]]
        }
    } else {
        # Standard matrix format
        expression_matrix <- as.matrix(read.table(data_dir, header = TRUE, row.names = 1))
    }

    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
        counts = expression_matrix,
        min.cells = min_cells,
        project = sample_id
    )

    # Add QC metrics
    seurat_obj[["percent.mt"]]   <- PercentageFeatureSet(seurat_obj, pattern = mito_grep)
    seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = ribo_grep)

    # Add sample metadata from samples file
    tryCatch({
        samples_meta <- read.csv("${samples_file}", header = TRUE, stringsAsFactors = FALSE)
        if ("condition" %in% colnames(samples_meta)) {
            cond <- samples_meta[samples_meta\$sample == sample_id, "condition"]
            if (length(cond) > 0) {
                seurat_obj\$condition <- cond[1]
            }
        }
    }, error = function(e) {
        message("Could not read samples metadata: ", e\$message)
    })

    # Generate QC plots
    pdf("${meta.id}_vlnplot_QC_prefilt.pdf", width = 12, height = 6)
    print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4))
    dev.off()

    pdf("${meta.id}_geneplot_numi_vs_pctmit_ngene.pdf", width = 12, height = 5)
    plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(plot1 + plot2)
    dev.off()

    # Save
    saveRDS(seurat_obj, file = "${meta.id}_seurat_pre-qc.rds")
    message("Seurat QC complete for ${meta.id}")
    """
}

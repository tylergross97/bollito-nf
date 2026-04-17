process SEURAT_MERGE {
    tag 'merged'
    label 'process_high_memory'

    conda 'conda-forge::r-seurat=5.1.0'

    input:
    path rds_files     // collected list of per-sample RDS files
    val sample_ids     // collected list of sample IDs
    val meta           // meta map passed from workflow

    output:
    tuple val(meta), path("merged_seurat_post-qc.rds"), emit: rds

    script:
    """
    #!/usr/bin/env Rscript

    suppressMessages(library("Seurat"))

    message("=== SEURAT MERGE ===")

    random_seed <- ${params.random_seed}
    set.seed(random_seed)
    write_table <- ${params.write_table ? 'TRUE' : 'FALSE'}

    # Load all sample objects
    rds_paths <- strsplit("${rds_files.join(',')}", ",")[[1]]
    sample_names <- strsplit("${sample_ids.join(',')}", ",")[[1]]

    objects <- list()
    for (i in seq_along(rds_paths)) {
        obj <- readRDS(trimws(rds_paths[i]))
        obj\$orig.sample <- sample_names[i]
        objects[[sample_names[i]]] <- obj
    }

    # Merge
    if (length(objects) > 1) {
        merged <- merge(objects[[1]], y = objects[2:length(objects)],
                        add.cell.ids = sample_names,
                        project = "merged")
    } else {
        merged <- objects[[1]]
    }

    message(paste0("Merged object: ", ncol(merged), " cells from ", length(objects), " samples"))

    saveRDS(merged, file = "merged_seurat_post-qc.rds")
    message("Merge complete")
    """
}

#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bollito-nf: single-cell RNA-seq analysis pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Migrated from the Snakemake bollito pipeline (CNIO Bioinformatics Unit).

    Supports:
      - FASTQ input (STARsolo alignment + QC) or pre-computed count matrices
      - 10x Chromium, Drop-seq, and custom scRNA-seq technologies
      - Seurat-based QC, normalization (standard / SCT), clustering, DEGs
      - Optional: merge, integration, gene-set scoring, trajectory inference,
        VISION functional analysis, RNA velocity

    Usage:
      nextflow run main.nf --input samplesheet.csv --fasta genome.fa --gtf genes.gtf [options]

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// Resolve a path that may be a local file or a remote (HTTP/FTP) URL.
// checkIfExists only works for local paths, not URLs.
def resolveFile(String path) {
    path ==~ /^https?:\/\/.*|^ftp:\/\/.*/ ? file(path) : file(path, checkIfExists: true)
}

// ── Include all process modules ─────────────────────────────────────────────
include { STAR_GENOMEGENERATE         } from './modules/local/star_genomegenerate'
include { STARSOLO                    } from './modules/local/starsolo'
include { FASTQC                      } from './modules/local/fastqc'
include { MULTIQC                     } from './modules/local/multiqc'
include { GTF2BED                     } from './modules/local/gtf2bed'
include { RSEQC_BAMSTAT               } from './modules/local/rseqc'
include { RSEQC_INFER_EXPERIMENT      } from './modules/local/rseqc'
include { RSEQC_JUNCTION_ANNOTATION   } from './modules/local/rseqc'
include { RSEQC_JUNCTION_SATURATION   } from './modules/local/rseqc'
include { RSEQC_INNER_DISTANCE        } from './modules/local/rseqc'
include { RSEQC_READ_DISTRIBUTION     } from './modules/local/rseqc'
include { RSEQC_READ_DUPLICATION      } from './modules/local/rseqc'
include { RSEQC_READ_GC              } from './modules/local/rseqc'
include { STAR_TO_VELOCYTO            } from './modules/local/star_to_velocyto'
include { SEURAT_QC                   } from './modules/local/seurat_qc'
include { SEURAT_POSTQC               } from './modules/local/seurat_postqc'
include { SEURAT_FILTER               } from './modules/local/seurat_filter'
include { SEURAT_MERGE                } from './modules/local/seurat_merge'
include { SEURAT_NORMALIZATION        } from './modules/local/seurat_normalization'
include { SEURAT_INTEGRATION          } from './modules/local/seurat_integration'
include { SEURAT_FIND_CLUSTERS        } from './modules/local/seurat_find_clusters'
include { SEURAT_DEGS                 } from './modules/local/seurat_degs'
include { SEURAT_GS_SCORING           } from './modules/local/seurat_gs_scoring'
include { SLINGSHOT                   } from './modules/local/slingshot'
include { VISION_ANALYSIS             } from './modules/local/vision_analysis'
include { RNA_VELOCITY                } from './modules/local/rna_velocity'

// ────────────────────────────────────────────────────────────────────────────
//  WORKFLOW
// ────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Validate required params ────────────────────────────────────────
    if (!params.input) { error "Please provide a samplesheet with --input" }

    // ── Parse samplesheet ───────────────────────────────────────────────
    // CSV format:
    //   FASTQ:            sample,condition,fq1,fq2
    //   Matrix (10x):     sample,condition,matrix_dir
    //   Matrix (standard):sample,condition,matrix_file
    ch_input_raw = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [id: row.sample, condition: row.condition ?: row.sample]
            if (params.input_type == 'fastq') {
                def fq1 = resolveFile(row.fq1)
                def fq2 = resolveFile(row.fq2)
                return [meta, fq1, fq2]
            } else {
                def mat = resolveFile(row.matrix_dir ?: row.matrix_file)
                return [meta, mat]
            }
        }

    ch_multiqc_files = channel.empty()

    // Reference files
    ch_gtf = params.gtf
        ? channel.of(resolveFile(params.gtf)).collect()
        : channel.empty()

    ch_samples_tsv = params.samples_tsv
        ? channel.fromPath(params.samples_tsv, checkIfExists: true).collect()
        : channel.fromPath(params.input, checkIfExists: true).collect()

    // Declare channels used across branches
    ch_seurat_input   = channel.empty()
    ch_velocyto_dirs  = channel.empty()

    // ====================================================================
    //  BRANCH A: FASTQ input → alignment + QC
    // ====================================================================
    if (params.input_type == 'fastq') {

        if (!params.fasta) { error "Please provide a reference genome with --fasta" }
        if (!params.gtf)   { error "Please provide a GTF annotation with --gtf" }

        ch_fasta = channel.of(resolveFile(params.fasta)).collect()

        // STAR genome index
        if (params.star_index) {
            ch_star_index = channel.fromPath(params.star_index, checkIfExists: true).collect()
        } else {
            STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
            ch_star_index = STAR_GENOMEGENERATE.out.index
        }

        // Whitelist
        ch_whitelist = params.whitelist
            ? channel.of(resolveFile(params.whitelist)).collect()
            : channel.value(file('NO_FILE'))

        // Group FASTQ by sample (merge multi-unit/lane rows)
        ch_reads_grouped = ch_input_raw
            .map { meta, fq1, fq2 -> [meta.id, meta, fq1, fq2] }
            .groupTuple(by: 0)
            .map { sample_id, metas, fq1s, fq2s ->
                def meta = metas[0]
                def all_reads = (fq1s + fq2s).flatten()
                [meta, all_reads]
            }

        // FastQC
        FASTQC(ch_reads_grouped)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { meta, zip -> zip })

        // STARsolo alignment
        STARSOLO(ch_reads_grouped, ch_star_index, ch_gtf, ch_whitelist)

        // RSeQC QC (optional)
        if (params.run_rseqc) {
            GTF2BED(ch_gtf)
            ch_bed = GTF2BED.out.bed.collect()

            RSEQC_BAMSTAT(STARSOLO.out.bam)
            RSEQC_INFER_EXPERIMENT(STARSOLO.out.bam, ch_bed)
            RSEQC_JUNCTION_ANNOTATION(STARSOLO.out.bam, ch_bed)
            RSEQC_JUNCTION_SATURATION(STARSOLO.out.bam, ch_bed)
            RSEQC_INNER_DISTANCE(STARSOLO.out.bam, ch_bed)
            RSEQC_READ_DISTRIBUTION(STARSOLO.out.bam, ch_bed)
            RSEQC_READ_DUPLICATION(STARSOLO.out.bam)
            RSEQC_READ_GC(STARSOLO.out.bam)

            ch_multiqc_files = ch_multiqc_files
                .mix(RSEQC_BAMSTAT.out.stats.map { meta, f -> f })
                .mix(RSEQC_INFER_EXPERIMENT.out.txt.map { meta, f -> f })
                .mix(RSEQC_JUNCTION_ANNOTATION.out.log.map { meta, f -> f })
                .mix(RSEQC_READ_DISTRIBUTION.out.txt.map { meta, f -> f })
        }

        ch_multiqc_files = ch_multiqc_files
            .mix(STARSOLO.out.log.map { meta, f -> f })

        // MultiQC
        ch_multiqc_config = params.multiqc_config
            ? channel.fromPath(params.multiqc_config, checkIfExists: true)
            : channel.value(file('NO_FILE'))

        MULTIQC(ch_multiqc_files.collect(), ch_multiqc_config)

        // Seurat input = Solo.out dirs
        ch_seurat_input = STARSOLO.out.solo_out

        // Velocyto matrices
        if (params.run_velocyto) {
            STAR_TO_VELOCYTO(STARSOLO.out.solo_out)
            ch_velocyto_dirs = STAR_TO_VELOCYTO.out.velocyto_dir
        }

    } else {
        // ================================================================
        //  BRANCH B: Matrix input → skip alignment
        // ================================================================
        ch_seurat_input = ch_input_raw
    }

    // ====================================================================
    //  SEURAT ANALYSIS PIPELINE (common for both input types)
    // ====================================================================

    // Step 1: QC
    SEURAT_QC(ch_seurat_input, ch_samples_tsv)

    // Step 2: Post-QC filtering
    SEURAT_POSTQC(SEURAT_QC.out.rds)

    // Step 2.1: Optional gene filtering
    ch_postqc_rds = SEURAT_POSTQC.out.rds

    if (params.filter_genes) {
        SEURAT_FILTER(SEURAT_POSTQC.out.rds)
        ch_postqc_rds = SEURAT_FILTER.out.rds
    }

    // Step 2.5: Merge (optional, multi-sample only)
    ch_merged_rds = channel.empty()

    if (params.run_merge) {
        ch_merge_rds_collected = ch_postqc_rds
            .map { meta, rds -> rds }
            .collect()
        ch_merge_ids_collected = ch_postqc_rds
            .map { meta, rds -> meta.id }
            .collect()

        SEURAT_MERGE(ch_merge_rds_collected, ch_merge_ids_collected, [id: 'merged', condition: 'merged'])
        ch_merged_rds = SEURAT_MERGE.out.rds
    }

    // Step 3: Normalization
    // Normalize individual samples (always) plus merged (if enabled)
    ch_to_normalize = params.single_samples
        ? ch_postqc_rds.mix(ch_merged_rds)
        : (params.run_merge ? ch_merged_rds : ch_postqc_rds)

    // Step 3B: Integration (optional, multi-sample only)
    ch_integration_rds = channel.empty()

    if (params.run_integration) {
        ch_integration_rds_collected = ch_postqc_rds
            .map { meta, rds -> rds }
            .collect()
        ch_integration_ids_collected = ch_postqc_rds
            .map { meta, rds -> meta.id }
            .collect()

        SEURAT_INTEGRATION(ch_integration_rds_collected, ch_integration_ids_collected, [id: 'integrated', condition: 'integrated'])
        ch_integration_rds = SEURAT_INTEGRATION.out.rds
    }

    SEURAT_NORMALIZATION(ch_to_normalize)

    // Step 4: Clustering
    // Cluster both normalized single/merged objects and integrated objects
    ch_to_cluster = SEURAT_NORMALIZATION.out.rds.mix(ch_integration_rds)

    SEURAT_FIND_CLUSTERS(ch_to_cluster)

    // Step 5: DEGs (optional)
    if (params.run_degs) {
        SEURAT_DEGS(SEURAT_FIND_CLUSTERS.out.rds)
    }

    // Step 6: Gene set scoring (optional)
    if (params.run_gs && params.geneset_collection) {
        ch_gs_gmt = channel.fromPath(params.geneset_collection, checkIfExists: true).collect()
        SEURAT_GS_SCORING(SEURAT_FIND_CLUSTERS.out.rds, ch_gs_gmt)
    }

    // Step 7: Trajectory inference (optional)
    if (params.run_slingshot) {
        SLINGSHOT(SEURAT_FIND_CLUSTERS.out.rds)
    }

    // Step 8: VISION functional analysis (optional)
    if (params.run_vision && (params.vision_geneset || params.geneset_collection)) {
        ch_vision_gmt = channel.fromPath(
            params.vision_geneset ?: params.geneset_collection, checkIfExists: true
        ).collect()
        VISION_ANALYSIS(SEURAT_FIND_CLUSTERS.out.rds, ch_vision_gmt)
    }

    // Step 9: RNA velocity (optional, FASTQ mode only)
    if (params.input_type == 'fastq' && params.run_velocyto) {
        // Join clustered Seurat objects with their velocyto matrices by meta.id
        ch_clustered_for_velocity = SEURAT_FIND_CLUSTERS.out.rds
            .filter { meta, rds -> meta.id != 'merged' && meta.id != 'integrated' }

        ch_velocity_joined = ch_clustered_for_velocity
            .map { meta, rds -> [meta.id, meta, rds] }
            .join(
                ch_velocyto_dirs.map { meta, vdir -> [meta.id, vdir] },
                by: 0
            )
            .map { sample_id, meta, rds, vdir -> [meta, rds, vdir] }

        RNA_VELOCITY(
            ch_velocity_joined.map { meta, rds, vdir -> [meta, rds] },
            ch_velocity_joined.map { meta, rds, vdir -> [meta, vdir] }
        )
    }
}

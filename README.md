# bollito-nf

A Nextflow DSL2 port of [bollito](https://gitlab.com/bu_cnio/bollito), the single-cell RNA-seq analysis pipeline from the CNIO Bioinformatics Unit.

## What Seqera AI did

Seqera AI converted the original Snakemake-based bollito pipeline into a fully modular Nextflow DSL2 pipeline. The conversion includes:

- **`main.nf`** — complete workflow definition with channel-based dataflow connecting all analysis steps
- **`nextflow.config`** — all pipeline parameters, resource labels, execution profiles (Docker, Singularity, Wave), and built-in timeline/trace reporting
- **`modules/local/`** — 20+ individual process modules, one per tool:
  - `starsolo.nf`, `star_genomegenerate.nf`, `star_to_velocyto.nf` — STAR alignment and index building
  - `fastqc.nf`, `multiqc.nf`, `rseqc.nf`, `gtf2bed.nf` — read QC and reporting
  - `seurat_qc.nf`, `seurat_postqc.nf`, `seurat_filter.nf` — cell/gene quality filtering
  - `seurat_merge.nf`, `seurat_integration.nf` — multi-sample merging and batch correction
  - `seurat_normalization.nf`, `seurat_find_clusters.nf`, `seurat_degs.nf` — normalization, clustering, differential expression
  - `seurat_gs_scoring.nf`, `slingshot.nf`, `vision_analysis.nf`, `rna_velocity.nf` — gene-set scoring, trajectory inference, VISION, RNA velocity
- **`assets/`** — example samplesheets for FASTQ and matrix inputs
- **`conf/test.config`** — lightweight test profile using publicly available reference data

Wave container support is enabled by default so no manual container management is needed.

## Pipeline overview

```
FASTQ input                       Matrix input
     │                                 │
     ▼                                 │
STAR_GENOMEGENERATE (optional)         │
     │                                 │
     ▼                                 │
STARsolo ──► FastQC ──► MultiQC        │
     │                                 │
     ▼                                 │
RSeQC (optional)                       │
     │                                 ▼
     └──────────────────► SEURAT_QC ──► SEURAT_POSTQC
                                │
                        SEURAT_FILTER (optional)
                                │
                     SEURAT_MERGE / SEURAT_INTEGRATION (optional)
                                │
                        SEURAT_NORMALIZATION
                                │
                        SEURAT_FIND_CLUSTERS
                                │
               ┌────────────────┼────────────────┐
               ▼                ▼                ▼
          SEURAT_DEGS    SEURAT_GS_SCORING   SLINGSHOT
               │                                 │
               ▼                                 ▼
        VISION_ANALYSIS                    RNA_VELOCITY
```

## Quick start

```bash
nextflow run main.nf \
  --input assets/samplesheet.csv \
  --fasta /path/to/genome.fa \
  --gtf   /path/to/annotation.gtf \
  -profile docker
```

For a minimal smoke-test using public data:

```bash
nextflow run main.nf -profile test,docker
```

## Input samplesheet

| Column | Description |
|---|---|
| `sample` | Unique sample name |
| `condition` | Experimental condition / group |
| `fq1` | Path or URL to R1 FASTQ (FASTQ mode) |
| `fq2` | Path or URL to R2 FASTQ (FASTQ mode) |
| `matrix_dir` | Path to 10x-format directory (matrix mode) |
| `matrix_file` | Path to standard matrix file (matrix mode) |

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_type` | `fastq` | `fastq` or `matrix` |
| `--technology` | `10x` | `10x`, `drop-seq`, `custom`, `standard` |
| `--technology_version` | `v3` | 10x chemistry version |
| `--norm_type` | `SCT` | `SCT` or `standard` |
| `--resolutions` | `0.2,0.4,0.8,1.2,1.6` | Seurat clustering resolutions |
| `--run_merge` | `true` | Merge multiple samples |
| `--run_integration` | `true` | Seurat integration / batch correction |
| `--run_degs` | `true` | Differential expression |
| `--run_slingshot` | `true` | Trajectory inference |
| `--run_vision` | `true` | VISION functional analysis |
| `--run_velocyto` | `true` | RNA velocity (FASTQ mode only) |

## Execution profiles

| Profile | Description |
|---|---|
| `docker` | Run with Docker |
| `singularity` | Run with Singularity |
| `wave` | Seqera Wave (auto-builds containers from conda; default) |
| `test` | Minimal test run with public data |

## Credits

Original pipeline: [CNIO Bioinformatics Unit — bollito](https://gitlab.com/bu_cnio/bollito)  
Nextflow conversion: Seqera AI

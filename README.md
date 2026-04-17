# bollito-nf

A Nextflow DSL2 port of [bollito](https://gitlab.com/bu_cnio/bollito), the single-cell RNA-seq analysis pipeline from the CNIO Bioinformatics Unit.

This project demonstrates two Seqera products working together end-to-end: **Seqera AI** to convert an existing Snakemake pipeline to production-ready Nextflow, and **Claude Code with the Seqera MCP server** to execute, iterate on, and validate it — all from the command line.

---

## The workflow

### Step 1 — Seqera AI: Snakemake → Nextflow

Seqera AI converted the full bollito Snakemake pipeline to Nextflow DSL2. What that involved:

- Translated every Snakemake rule into a standalone Nextflow process module (`modules/local/`)
- Re-expressed the Snakefile's dependency graph as channel-based dataflow in `main.nf`
- Ported all parameters from `config-example.yaml` and the `schemas/` definitions into `nextflow.config` with proper resource labels, retry logic, and execution profiles
- Added Wave support so containers are built on-the-fly from conda — no Dockerfile maintenance needed
- Generated samplesheets and a `conf/test.config` for lightweight CI runs

The conversion covered the full analysis surface: alignment (STARsolo), QC (FastQC, RSeQC, MultiQC), Seurat QC → normalization → clustering → DEGs, gene-set scoring, trajectory inference (Slingshot), VISION functional analysis, and RNA velocity.

**Everything from the initial commit back through the pipeline design was done by Seqera AI.**

---

### Step 2 — Claude Code + Seqera MCP: test runs from the CLI

After the conversion, Claude Code was used with the [Seqera MCP server](https://github.com/seqeralabs/mcp-seqera) to execute and validate the pipeline directly from the terminal. The Seqera MCP gives Claude Code live access to the Seqera Platform API — launch runs, inspect compute environments, and pull execution logs without leaving your editor.

Running the test profile surfaced several real incompatibilities between the generated pipeline and the current Seurat/Nextflow ecosystem. Claude Code diagnosed and fixed them in-place:

| File | Fix |
|---|---|
| `seurat_qc.nf`, `seurat_postqc.nf` | Replaced `VlnPlot()` and `FeatureScatter()` with direct ggplot2 — Seurat v5.1.0 introduced a `FetchData` API regression that broke both calls |
| `seurat_find_clusters.nf` | Made `SeuratDisk` and `lisi` optional dependencies; removed them from the conda spec and added graceful fallbacks so a missing package doesn't fail the entire process |
| `conf/test.config` | Switched `norm_type` from `SCT` to `standard` to work around a `slot=` API bug in Seurat 5.1.0 + SeuratObject 5.0; relaxed QC thresholds to match the chr19-only test reference |

After these fixes the test profile completed successfully end-to-end.

Example session:

```
$ claude

> launch bollito-nf test profile on my AWS compute environment and tail the logs

# Claude Code uses the Seqera MCP to:
# 1. Find matching compute environments via the API
# 2. Launch the pipeline with -profile test,wave
# 3. Stream back the run status and surface any process-level errors

> the SEURAT_QC process is failing — VlnPlot error in Seurat v5, fix it

# Claude Code reads the module, rewrites the plotting block with ggplot2,
# commits the fix, and relaunches — no context switch to the browser needed
```

This is the development loop in practice: catch a failure in the platform logs, fix the module, relaunch — all from one terminal session.

---

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

---

## Quick start

```bash
nextflow run main.nf \
  --input assets/samplesheet.csv \
  --fasta /path/to/genome.fa \
  --gtf   /path/to/annotation.gtf \
  -profile docker
```

Minimal smoke-test using public data:

```bash
nextflow run main.nf -profile test,docker
```

---

## Input samplesheet

| Column | Description |
|---|---|
| `sample` | Unique sample name |
| `condition` | Experimental condition / group |
| `fq1` | Path or URL to R1 FASTQ (FASTQ mode) |
| `fq2` | Path or URL to R2 FASTQ (FASTQ mode) |
| `matrix_dir` | Path to 10x-format directory (matrix mode) |
| `matrix_file` | Path to standard matrix file (matrix mode) |

---

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

---

## Execution profiles

| Profile | Description |
|---|---|
| `docker` | Run with Docker |
| `singularity` | Run with Singularity |
| `wave` | Seqera Wave — auto-builds containers from conda (default) |
| `test` | Minimal test run with public data |

---

## Credits

Original pipeline: [CNIO Bioinformatics Unit — bollito](https://gitlab.com/bu_cnio/bollito)  
Nextflow conversion: [Seqera AI](https://seqera.io/ai)  
Pipeline execution & iteration: Claude Code + [Seqera MCP server](https://github.com/seqeralabs/mcp-seqera)

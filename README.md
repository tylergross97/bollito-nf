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

After the conversion, Claude Code was used with the [Seqera MCP server](https://github.com/seqeralabs/mcp-seqera) to run and iterate on the pipeline directly from the terminal — no browser required.

The Seqera MCP server gives Claude Code direct access to the Seqera Platform API, so you can launch runs, inspect compute environments, manage datasets, and pull execution logs without leaving your editor.

Example session:

```
$ claude

> launch bollito-nf test profile on my AWS compute environment and tail the logs

# Claude Code uses the Seqera MCP to:
# 1. Find matching compute environments via the API
# 2. Launch the pipeline with -profile test
# 3. Stream back the run status and any process-level errors
```

This makes the development loop much tighter: catch a process failure, fix the module, relaunch — without switching context to the Seqera Platform UI.

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

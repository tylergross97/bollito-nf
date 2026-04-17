process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.zip"), emit: zip
    tuple val(meta), path("*.html"), emit: html

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

process STAR_GENOMEGENERATE {
    tag 'genome'
    label 'process_high'

    conda 'bioconda::star=2.7.11b'

    input:
    path fasta
    path gtf

    output:
    path 'star_index', emit: index

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --runThreadN ${task.cpus} \\
        ${args}
    """
}

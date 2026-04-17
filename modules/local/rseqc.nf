process RSEQC_BAMSTAT {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::rseqc=5.0.4'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats.txt"), emit: stats

    script:
    """
    bam_stat.py -i ${bam} > ${meta.id}.stats.txt
    """
}

process RSEQC_INFER_EXPERIMENT {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::rseqc=5.0.4'

    input:
    tuple val(meta), path(bam)
    path bed

    output:
    tuple val(meta), path("*.infer_experiment.txt"), emit: txt

    script:
    """
    infer_experiment.py -r ${bed} -i ${bam} > ${meta.id}.infer_experiment.txt
    """
}

process RSEQC_JUNCTION_ANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::rseqc=5.0.4'

    input:
    tuple val(meta), path(bam)
    path bed

    output:
    tuple val(meta), path("*.junction.bed"), emit: bed_out
    tuple val(meta), path("*.log"),          emit: log

    script:
    """
    junction_annotation.py -q 255 -i ${bam} -r ${bed} -o ${meta.id}.junctionanno 2>&1 | tee ${meta.id}.junctionanno.log
    """
}

process RSEQC_JUNCTION_SATURATION {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::rseqc=5.0.4'

    input:
    tuple val(meta), path(bam)
    path bed

    output:
    tuple val(meta), path("*.pdf"), emit: pdf

    script:
    """
    junction_saturation.py -q 255 -i ${bam} -r ${bed} -o ${meta.id}.junctionsat
    """
}

process RSEQC_INNER_DISTANCE {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::rseqc=5.0.4'

    input:
    tuple val(meta), path(bam)
    path bed

    output:
    tuple val(meta), path("*.inner_distance.txt"), emit: txt

    script:
    """
    inner_distance.py -r ${bed} -i ${bam} -o ${meta.id}.inner_distance_freq
    """
}

process RSEQC_READ_DISTRIBUTION {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::rseqc=5.0.4'

    input:
    tuple val(meta), path(bam)
    path bed

    output:
    tuple val(meta), path("*.readdistribution.txt"), emit: txt

    script:
    """
    read_distribution.py -r ${bed} -i ${bam} > ${meta.id}.readdistribution.txt
    """
}

process RSEQC_READ_DUPLICATION {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::rseqc=5.0.4'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf

    script:
    """
    read_duplication.py -i ${bam} -o ${meta.id}.readdup
    """
}

process RSEQC_READ_GC {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::rseqc=5.0.4'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf

    script:
    """
    read_GC.py -i ${bam} -o ${meta.id}.readgc
    """
}

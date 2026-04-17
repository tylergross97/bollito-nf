process STARSOLO {
    tag "${meta.id}"
    label 'process_high'

    conda 'bioconda::star=2.7.11b bioconda::samtools=1.21'

    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    path whitelist

    output:
    tuple val(meta), path("${meta.id}/Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}/Solo.out"),                      emit: solo_out
    tuple val(meta), path("${meta.id}/Log.final.out"),                 emit: log

    script:
    // Build technology-specific STAR params
    def star_extra = params.star_extra_args
    if (!star_extra) {
        if (params.technology == '10x') {
            def umi_len = params.technology_version == 'v1' ? 10 : (params.technology_version == 'v2' ? 10 : 12)
            def cb_len  = params.technology_version == 'v1' ? 14 : 16
            star_extra = "--soloType Droplet --soloFeatures Gene Velocyto --outFilterMultimapNmax 50 --winAnchorMultimapNmax 50 --alignEndsType EndToEnd --outReadsUnmapped Fastx --soloUMIlen ${umi_len} --soloCBlen ${cb_len} --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM"
        } else if (params.technology == 'drop-seq') {
            star_extra = "--soloType Droplet --soloFeatures Gene Velocyto --outFilterMultimapNmax 50 --winAnchorMultimapNmax 50 --alignEndsType EndToEnd --outReadsUnmapped Fastx --soloUMIlen 8 --soloCBlen 12 --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM"
        } else {
            star_extra = "--soloType Droplet --soloFeatures Gene Velocyto --outFilterMultimapNmax 50 --winAnchorMultimapNmax 50 --alignEndsType EndToEnd --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM"
        }
    }

    // Separate R1 (biological) and R2 (barcode) reads
    // STARsolo expects barcode read as the second file
    def read_list = reads.collect { it.name }
    def r1_files = read_list.findAll { it.contains('.r1.') || it.contains('_R1') || it.contains('_1.') }.sort().join(',')
    def r2_files = read_list.findAll { it.contains('.r2.') || it.contains('_R2') || it.contains('_2.') }.sort().join(',')

    // STARsolo expects: --readFilesIn cDNA_read barcode_read
    // For 10x: R1 = barcode+UMI, R2 = cDNA insert
    def read_files_in = "${r2_files} ${r1_files}"

    def whitelist_opt = whitelist.name != 'NO_FILE' ? "--soloCBwhitelist ${whitelist}" : "--soloCBwhitelist None"
    def gz_opt = reads[0].name.endsWith('.gz') ? "--readFilesCommand zcat" : ""

    """
    mkdir -p ${meta.id}
    STAR \\
        --genomeDir ${index} \\
        --sjdbGTFfile ${gtf} \\
        --readFilesIn ${read_files_in} \\
        ${whitelist_opt} \\
        ${gz_opt} \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${meta.id}/ \\
        ${star_extra}

    samtools index ${meta.id}/Aligned.sortedByCoord.out.bam

    # Compress STARsolo matrix files so Read10X can find barcodes.tsv.gz
    find ${meta.id}/Solo.out -name "*.tsv" -o -name "*.mtx" | xargs gzip --force
    """
}

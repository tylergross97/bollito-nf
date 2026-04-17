process MULTIQC {
    tag 'multiqc'
    label 'process_single'

    conda 'bioconda::multiqc=1.25.2'

    input:
    path('multiqc_input/*')
    path multiqc_config

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data',        emit: data

    script:
    def config_opt = multiqc_config.name != 'NO_FILE' ? "--config ${multiqc_config}" : ""
    """
    multiqc ${config_opt} --force multiqc_input/
    """
}

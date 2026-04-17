process GTF2BED {
    tag 'gtf2bed'
    label 'process_single'

    conda 'conda-forge::python=3.12 bioconda::gffutils=0.13'

    input:
    path gtf

    output:
    path '*.bed', emit: bed

    script:
    """
    #!/usr/bin/env python3
    import gffutils

    db = gffutils.create_db(
        '${gtf}',
        dbfn='annotation.db',
        force=True,
        keep_order=True,
        merge_strategy='merge',
        sort_attribute_values=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True
    )

    with open('annotation.bed', 'w') as outfileobj:
        for tx in db.features_of_type('transcript', order_by='start'):
            bed = [s.strip() for s in db.bed12(tx).split('\\t')]
            bed[3] = tx.id
            outfileobj.write('{}\\n'.format('\\t'.join(bed)))
    """
}

process STAR_TO_VELOCYTO {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(solo_out)

    output:
    tuple val(meta), path("velocyto_out"), emit: velocyto_dir

    script:
    """
    mkdir -p velocyto_out

    INPUT_DIR="${solo_out}/Velocyto/raw"

    head -n 2 \${INPUT_DIR}/matrix.mtx > \${INPUT_DIR}/mtx_header.txt
    head -n 3 \${INPUT_DIR}/matrix.mtx | tail -n 1 > \${INPUT_DIR}/mtx_summary.txt

    tail -n +4 \${INPUT_DIR}/matrix.mtx | cut -d " " -f 1-3 > \${INPUT_DIR}/ms.txt
    tail -n +4 \${INPUT_DIR}/matrix.mtx | cut -d " " -f 1-2,4 > \${INPUT_DIR}/mu.txt
    tail -n +4 \${INPUT_DIR}/matrix.mtx | cut -d " " -f 1-2,5 > \${INPUT_DIR}/ma.txt

    mkdir -p velocyto_out/spliced velocyto_out/unspliced velocyto_out/ambiguous

    cat \${INPUT_DIR}/mtx_header.txt \${INPUT_DIR}/mtx_summary.txt \${INPUT_DIR}/ms.txt > velocyto_out/spliced/matrix.mtx
    cat \${INPUT_DIR}/mtx_header.txt \${INPUT_DIR}/mtx_summary.txt \${INPUT_DIR}/mu.txt > velocyto_out/unspliced/matrix.mtx
    cat \${INPUT_DIR}/mtx_header.txt \${INPUT_DIR}/mtx_summary.txt \${INPUT_DIR}/ma.txt > velocyto_out/ambiguous/matrix.mtx

    cp \${INPUT_DIR}/features.tsv velocyto_out/genes.tsv 2>/dev/null || true
    cp \${INPUT_DIR}/genes.tsv velocyto_out/genes.tsv 2>/dev/null || true
    cp velocyto_out/genes.tsv velocyto_out/spliced/
    cp velocyto_out/genes.tsv velocyto_out/unspliced/
    cp velocyto_out/genes.tsv velocyto_out/ambiguous/

    cp \${INPUT_DIR}/barcodes.tsv velocyto_out/
    cp velocyto_out/barcodes.tsv velocyto_out/spliced/
    cp velocyto_out/barcodes.tsv velocyto_out/unspliced/
    cp velocyto_out/barcodes.tsv velocyto_out/ambiguous/

    rm -f \${INPUT_DIR}/*.txt
    """
}

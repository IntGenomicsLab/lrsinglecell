process DROP_BARCODE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.10.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"

    input:
    // Staged under input/ so that the filtered files, written to the task directory, keep
    // the same names; the output globs do not reach into input/
    tuple val(meta), path(mtx, stageAs: 'input/*'), path(barcodes, stageAs: 'input/*'), path(features, stageAs: 'input/*')
    val barcode // str: the barcode to remove, e.g. flexiplex's '-' no-match placeholder

    output:
    tuple val(meta), path("*.{mtx,mtx.gz}")         , emit: mtx
    tuple val(meta), path("*barcodes.{tsv,tsv.gz}") , emit: barcodes
    tuple val(meta), path("*features.{tsv,tsv.gz}") , emit: features
    path "versions.yml"                             , emit: versions_drop_barcode, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    drop_barcode.py mex \\
        ${args} \\
        --barcode '${barcode}' \\
        --matrix ${mtx} \\
        --barcodes ${barcodes} \\
        --features ${features} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch \$(basename ${mtx}) \$(basename ${barcodes}) \$(basename ${features})

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process FLEXIPLEX_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flexiplex:1.02.7--py312h9b2995c_0':
        'biocontainers/flexiplex:1.02.7--py312h9b2995c_0' }"

    input:
    tuple val(meta), path(barcodes)
    path(whitelist)

    output:
    tuple val(meta), path("*known_barcodes.txt")        , emit: barcodes
    path "versions.yml"                                 , emit: versions_flexiplex_filter, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARN: the version is not taken from the binary. The bioconda flexiplex 1.02.7
    // build was made from the v1.02.7 tag before upstream bumped its version string,
    // so `flexiplex --help` still prints 'FLEXIPLEX 1.02.6.1'. Update this string
    // whenever the container is bumped.
    def VERSION = '1.02.7'
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}${meta.part ? "_part_${meta.part}" : ''}"
    """
    flexiplex-filter \\
        ${barcodes} \\
        --whitelist ${whitelist} \\
        --outfile ${prefix}_known_barcodes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: ${VERSION}
    END_VERSIONS
    """

    stub:
    def VERSION = '1.02.7' // WARN: see the script block
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_known_barcodes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: ${VERSION}
    END_VERSIONS
    """
}

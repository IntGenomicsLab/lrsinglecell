process SAMTOOLS_FILTER_SORT {
    tag "$meta.id"
    label 'process_medium'

    // Same samtools build and container as the nf-core samtools/sort module
    conda "bioconda::samtools=1.22.1 bioconda::htslib=1.22.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml"                   , emit: versions_samtools_filter_sort, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''  // samtools view: the filter
    def args2 = task.ext.args2 ?: ''  // samtools sort: the sort order
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    samtools view \\
        -u \\
        ${args} \\
        ${bam} \\
    | \\
    samtools sort \\
        ${args2} \\
        -T ${prefix} \\
        --threads ${task.cpus} \\
        -o ${prefix}.bam \\
        -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | sed '1!d;s/.* //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | sed '1!d;s/.* //')
    END_VERSIONS
    """
}

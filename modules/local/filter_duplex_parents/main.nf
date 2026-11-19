process FILTER_DUPLEX_PARENTS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8':
        'biocontainers/pigz:2.8' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.duplex_filtered.fastq.gz"), emit: reads
    tuple val(meta), path("${prefix}.duplex_filter.tsv")       , emit: stats
    tuple val(meta), env('mode')                               , emit: mode
    path "versions.yml"                                        , emit: versions_filter_duplex_parents, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // dorado duplex tags every FASTQ record with dx:i:1 (duplex consensus read),
    // dx:i:0 (simplex read without a partner) or dx:i:-1 (simplex parent of a
    // dx:i:1 read). The parents are dropped so each duplex molecule enters the
    // pipeline once. Tags may be tab- or space-separated.
    """
    out=${prefix}.duplex_filtered.fastq.gz

    # dorado writes dx on every record of a duplex-mode file, so the headers of
    # the first 1000 records tell whether this is one. pigz -dc | head ends with
    # SIGPIPE, which must not fail the task under pipefail.
    set +o pipefail
    pigz -dc "${reads}" | head -n 4000 | awk 'NR % 4 == 1' > probe_headers.txt || true
    set -o pipefail
    n_probe_dx=\$(grep -cE '(^|[[:space:]])dx:i:' probe_headers.txt || true)
    # Duplex consensus reads are named <template>;<complement>
    n_probe_pairs=\$(awk '\$1 ~ /;/ { n++ } END { print n + 0 }' probe_headers.txt)

    if [ "\${n_probe_dx:-0}" -eq 0 ] && [ "\$n_probe_pairs" -gt 0 ]; then
        # A duplex FASTQ without dx tags (dorado duplex --emit-fastq): the
        # parents cannot be told apart, so fail rather than silently keep them
        echo "ERROR: ${meta.id}: ${reads} contains dorado duplex consensus reads (read names with ';')" \\
            "but no dx tags, so the duplex parents cannot be identified." \\
            "'dorado duplex --emit-fastq' does not write dx. Convert dorado's duplex BAM instead," \\
            "keeping the tags: samtools fastq -T '*' calls.bam | gzip > reads.fastq.gz." \\
            "Or run without --drop_duplex_parents to keep all reads." >&2
        exit 1
    elif [ "\${n_probe_dx:-0}" -eq 0 ]; then
        # Not a dorado duplex FASTQ: pass it through without reading it in full
        mode=passthrough_no_dx
        n_in=NA; n_1=NA; n_0=NA; n_m1=NA; n_none=NA; dropped=0
        ln -s "${reads}" "\$out"
    else
        pigz -dc "${reads}" \\
            | awk -v counts=counts.txt '
                NR % 4 == 1 {
                    n_in++
                    drop = 0
                    if (\$0 ~ /(^|[[:space:]])dx:i:-1([[:space:]]|\$)/)     { drop = 1; n_m1++ }
                    else if (\$0 ~ /(^|[[:space:]])dx:i:1([[:space:]]|\$)/) { n_1++ }
                    else if (\$0 ~ /(^|[[:space:]])dx:i:0([[:space:]]|\$)/) { n_0++ }
                    else                                                    { n_none++ }
                }
                !drop
                END { printf "%d %d %d %d %d\\n", n_in, n_1, n_0, n_m1, n_none > counts }
            ' \\
            | pigz -p ${task.cpus} > "\$out"

        read n_in n_1 n_0 n_m1 n_none < counts.txt

        if [ "\$n_m1" -eq 0 ]; then
            mode=nothing_to_drop
            dropped=0
        elif [ "\$n_1" -eq 0 ]; then
            # Parents without any consensus read: the consensus reads were most
            # likely removed upstream, so dropping the parents would lose molecules
            mode=kept_all_no_consensus
            dropped=0
            echo "WARNING: ${meta.id}: \$n_m1 dx:i:-1 records but no dx:i:1 records in ${reads}." \\
                "Keeping all reads; dropping the duplex parents would lose these molecules." >&2
        else
            mode=filtered
            dropped=\$n_m1
        fi

        if [ "\$mode" != filtered ]; then
            # Keep the input untouched
            rm "\$out"
            ln -s "${reads}" "\$out"
        fi
    fi

    {
        printf 'metric\\tvalue\\n'
        printf 'sample\\t%s\\n'      "${meta.id}"
        printf 'mode\\t%s\\n'        "\$mode"
        printf 'records_in\\t%s\\n'  "\$n_in"
        printf 'dx_1\\t%s\\n'        "\$n_1"
        printf 'dx_0\\t%s\\n'        "\$n_0"
        printf 'dx_minus1\\t%s\\n'   "\$n_m1"
        printf 'no_dx_tag\\t%s\\n'   "\$n_none"
        printf 'dropped\\t%s\\n'     "\$dropped"
    } > ${prefix}.duplex_filter.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mode=stub
    echo '' | gzip > ${prefix}.duplex_filtered.fastq.gz
    printf 'metric\\tvalue\\n' > ${prefix}.duplex_filter.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/^pigz //' )
    END_VERSIONS
    """
}

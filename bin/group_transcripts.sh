#!/usr/bin/env bash

# Originally written by Austyn Trull (nf-core/scnanoseq) and released under the MIT license.
# Rewritten for IntGenomicsLab/lrsinglecell.
#
# Group the sequences of a transcriptome FASTA by the chromosome their transcript lies on,
# so a transcriptome BAM can be split into one chunk per chromosome.
#
# For every chromosome it writes `<chr>.transcripts.txt`, one BAM reference name per line, in
# FASTA order. The reference name is the first whitespace-delimited token of the FASTA header,
# which is what minimap2 (and so `read.reference_name` in split_bam.py) uses. Sequences whose
# transcript is not found in the GTF go to `chr_Unknown.transcripts.txt`.
#
# The transcript id is the part of the reference name before the delimiter (-d). Without -d,
# the delimiter is taken from the first header that contains `|` or `;` (whichever occurs most),
# and a reference name without one is used whole. The id is looked up in a transcript_id ->
# seqname map parsed from the `transcript_id "..."` attribute of GTF column 9, first as is and
# then with a trailing `.N` version removed on either side, so versioned and unversioned FASTA
# and GTF ids match each other.
#
# Usage: group_transcripts.sh -f transcripts.fa -g annotation.gtf [-d '|']

set -euo pipefail

########################
### INPUT PARAMETERS ###
########################

IN_FASTA=""
IN_GTF=""
DELIMITER=""

usage() {
    echo "Usage: $(basename "$0") -f <transcripts.fa> -g <annotation.gtf> [-d <delimiter>]" >&2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f) IN_FASTA="${2:-}"; shift;;
        -g) IN_GTF="${2:-}"; shift;;
        -d) DELIMITER="${2:-}"; shift;;
        -h|--help) usage; exit 0;;
        *) echo "ERROR: unknown argument '$1'" >&2; usage; exit 1;;
    esac
    shift
done

if [[ -z "$IN_FASTA" || -z "$IN_GTF" ]]; then
    usage
    exit 1
fi

for f in "$IN_FASTA" "$IN_GTF"; do
    if [[ ! -r "$f" ]]; then
        echo "ERROR: cannot read '$f'" >&2
        exit 1
    fi
done

#########################
### GROUP TRANSCRIPTS ###
#########################

# Pass 1 (the GTF) builds the map; pass 2 (the FASTA) looks every header up in it.
# Output files are opened once and appended to, so the run is linear in input size.

awk -F '\t' -v delim="$DELIMITER" '
    function strip_version(id) {
        sub(/\.[0-9]+$/, "", id)
        return id
    }

    # Keep the lexicographically first seqname if a transcript id occurs on several
    # (e.g. PAR copies), so the result does not depend on GTF order.
    function remember(map, id, seq) {
        if (!(id in map) || seq < map[id]) map[id] = seq
    }

    # Anything outside a conservative character set would make an awkward file name
    function file_safe(name) {
        gsub(/[^A-Za-z0-9._-]/, "_", name)
        return name
    }

    FNR == NR {
        sub(/\r$/, "")
        if ($0 ~ /^#/ || NF < 9) next
        if (!match($9, /(^|[ \t;])transcript_id[ \t]+"[^"]*"/)) next
        id = substr($9, RSTART, RLENGTH)
        sub(/^[ \t;]*transcript_id[ \t]+"/, "", id)
        sub(/"$/, "", id)
        if (id == "") next
        seq = $1 ""
        remember(by_id, id, seq)
        remember(by_base, strip_version(id), seq)
        n_gtf++
        next
    }

    /^>/ {
        sub(/\r$/, "")
        ref = substr($0, 2)
        sub(/^[ \t]+/, "", ref)
        sub(/[ \t].*$/, "", ref)
        if (ref == "") next

        if (delim == "") {
            tmp = ref; n_pipe = gsub(/\|/, "", tmp)
            tmp = ref; n_semi = gsub(/;/, "", tmp)
            if (n_pipe > 0 && n_pipe >= n_semi) delim = "|"
            else if (n_semi > 0) delim = ";"
        }

        key = ref
        if (delim != "" && (p = index(ref, delim)) > 0) key = substr(ref, 1, p - 1)
        base = strip_version(key)

        if (key in by_id)         chr = by_id[key]
        else if (base in by_id)   chr = by_id[base]
        else if (key in by_base)  chr = by_base[key]
        else if (base in by_base) chr = by_base[base]
        else                      chr = "chr_Unknown"

        # Truncate each output once, then append. Files are closed in batches so a
        # GTF with thousands of scaffolds cannot run into the open-file limit.
        out = file_safe(chr) ".transcripts.txt"
        if (!(out in seen)) {
            seen[out] = 1
            printf "" > out
            close(out)
        }
        print ref >> out
        if (!(out in is_open)) {
            is_open[out] = 1
            if (++n_open >= 200) {
                for (o in is_open) close(o)
                split("", is_open)
                n_open = 0
            }
        }
        n_fasta++
        if (chr == "chr_Unknown") n_unknown++
    }

    END {
        if (n_gtf == 0) {
            print "ERROR: no transcript_id attribute found in the GTF" > "/dev/stderr"
            exit 1
        }
        if (n_fasta == 0) {
            print "ERROR: no sequence headers found in the FASTA" > "/dev/stderr"
            exit 1
        }
        if (n_unknown == n_fasta) {
            print "WARNING: none of the " n_fasta " FASTA sequences matched a GTF transcript_id; all are in chr_Unknown.transcripts.txt" > "/dev/stderr"
        } else if (n_unknown > 0) {
            print "WARNING: " n_unknown " of " n_fasta " FASTA sequences matched no GTF transcript_id and are in chr_Unknown.transcripts.txt" > "/dev/stderr"
        }
    }
' "$IN_GTF" "$IN_FASTA"

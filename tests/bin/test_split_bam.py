#!/usr/bin/env python3

""" Unit test for bin/split_bam.py: the reads it keeps and the counts it reports.

    Not run by nf-test (which only collects *.nf.test); .github/workflows/bin-tests.yml
    runs it on pull requests. Locally, run it in the module's container, which has pysam
    but no pytest:

        apptainer exec <mulled-v2-bb96c7354781ab52d8e69ccff89587598dc87fea image> \
            python -m unittest tests/bin/test_split_bam.py
"""

import os
import subprocess
import sys
import tempfile
import unittest

import pysam

pysam.set_verbosity(0)

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin", "split_bam.py")

# Reference names as minimap2 writes them: the first token of each FASTA header
REFERENCES = [
    "ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|",
    "ENST00000473358.1",
    "ENST00000361390.2",
]


class SplitBamTest(unittest.TestCase):

    def setUp(self):
        self.workdir = tempfile.mkdtemp(prefix="split_bam_", dir=os.environ.get("TMPDIR"))
        self.addCleanup(subprocess.run, ["rm", "-rf", self.workdir])
        self.bam = os.path.join(self.workdir, "in.bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": name, "LN": 2000} for name in REFERENCES]}
        with pysam.AlignmentFile(self.bam, "wb", header=header) as out:
            # two reads on each reference
            for tid in range(len(REFERENCES)):
                for i in range(2):
                    read = pysam.AlignedSegment(out.header)
                    read.query_name = f"r{tid}_{i}"
                    read.reference_id = tid
                    read.reference_start = 100 + i
                    read.mapping_quality = 60
                    read.cigarstring = "4M"
                    read.query_sequence = "ACGT"
                    read.query_qualities = pysam.qualitystring_to_array("IIII")
                    out.write(read)
        pysam.index(self.bam)

    def split(self, contigs):
        contig_file = os.path.join(self.workdir, "contigs.txt")
        with open(contig_file, "w") as f:
            f.write("".join(c + "\n" for c in contigs))
        out_bam = os.path.join(self.workdir, "out.split.bam")
        counts = os.path.join(self.workdir, "out.read_counts.tsv")
        subprocess.run([sys.executable, SCRIPT, "-i", self.bam, "-o", out_bam, "-c", contig_file,
                        "--counts", counts], check=True)
        with pysam.AlignmentFile(out_bam, "rb") as f:
            names = sorted(r.query_name for r in f)
        with open(counts) as f:
            n_input, n_written = (int(x) for x in f.read().split("\t"))
        return names, n_input, n_written

    def test_keeps_reads_on_listed_references(self):
        names, n_input, n_written = self.split([REFERENCES[0], REFERENCES[2]])
        self.assertEqual(names, ["r0_0", "r0_1", "r2_0", "r2_1"])
        self.assertEqual((n_input, n_written), (6, 4))

    def test_full_header_matches_nothing(self):
        # What the old group_transcripts.sh wrote for an Ensembl cDNA header: the whole
        # line, which is never a BAM reference name. The counts let the caller notice.
        names, n_input, n_written = self.split(["ENST00000473358.1 cdna chromosome:GRCh38:1:29554:31097:1"])
        self.assertEqual(names, [])
        self.assertEqual((n_input, n_written), (6, 0))

    def test_counts_are_optional(self):
        contig_file = os.path.join(self.workdir, "contigs.txt")
        with open(contig_file, "w") as f:
            f.write(REFERENCES[1] + "\n")
        out_bam = os.path.join(self.workdir, "out.split.bam")
        subprocess.run([sys.executable, SCRIPT, "-i", self.bam, "-o", out_bam, "-c", contig_file], check=True)
        with pysam.AlignmentFile(out_bam, "rb") as f:
            self.assertEqual(sorted(r.query_name for r in f), ["r1_0", "r1_1"])


if __name__ == "__main__":
    unittest.main()

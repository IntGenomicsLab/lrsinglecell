#!/usr/bin/env python3

""" Unit test for bin/tag_genes.py: gene tagging, the bam/gtf contig name check and the
    per-sample GS=none summary (--tag_summary).

    Not run by nf-test (which only collects *.nf.test); .github/workflows/bin-tests.yml
    runs it on pull requests. Locally, run it in the module's container, which has pysam
    but no pytest:

        apptainer exec <pysam:0.19.1--py310hff46b53_1 image> \
            python -m unittest tests/bin/test_tag_genes.py
"""

import os
import subprocess
import sys
import tempfile
import unittest

import pysam

pysam.set_verbosity(0)

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin", "tag_genes.py")

# contig, start, end (1-based inclusive), gene_id, gene_name
GENES = [
    ("chr1", 1001, 2000, "G1", "GeneOne"),
    ("chr1", 5001, 6000, "G2", "GeneTwo"),
    # overlaps G2 over 5501-6000, so a read on 5600-5800 is ambiguous
    ("chr1", 5501, 7000, "G3", "GeneThree"),
]

# name, contig, 0-based start, aligned length; the expected GS for the chr1 gtf
READS = [
    ("in_g1", "chr1", 1100, 200, "unique"),
    ("in_g2_g3", "chr1", 5600, 200, "ambiguous"),
    ("intergenic", "chr1", 3000, 200, "none"),
    ("on_chr2", "chr2", 1000, 200, "none"),
]


def write_gtf(path, genes):
    with open(path, "w") as fh:
        fh.write("#!genome-build toy\n")
        for contig, start, end, gene_id, gene_name in genes:
            attrs = f'gene_id "{gene_id}"; gene_name "{gene_name}";'
            fh.write(f"{contig}\ttoy\tgene\t{start}\t{end}\t.\t+\t.\t{attrs}\n")
            # a non-gene row the script has to skip
            fh.write(f"{contig}\ttoy\texon\t{start}\t{end}\t.\t+\t.\t{attrs} transcript_id \"T{gene_id}\";\n")


def renamed(genes, rename):
    return [(rename(contig), *rest) for contig, *rest in genes]


def build_bam(path, contigs=("chr1", "chr2")):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": contig, "LN": 100000} for contig in contigs],
    })
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for name, contig, start, length, _ in sorted(READS, key=lambda r: (contigs.index(r[1]), r[2])):
            r = pysam.AlignedSegment(header)
            r.query_name = name
            r.reference_name = contig
            r.reference_start = start
            r.mapping_quality = 60
            r.cigartuples = [(0, length)]
            r.query_sequence = "A" * length
            r.query_qualities = pysam.qualitystring_to_array("I" * length)
            out.write(r)
    pysam.index(path)


class TagGenesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.bam = os.path.join(cls.tmp.name, "in.bam")
        build_bam(cls.bam)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_tool(self, tag, genes, extra_args=()):
        d = self.tmp.name
        gtf = os.path.join(d, f"{tag}.gtf")
        write_gtf(gtf, genes)
        out = {k: os.path.join(d, f"{tag}.{k}") for k in ("bam", "tsv")}
        out["tag_summary"] = os.path.join(d, f"{tag}.gene_tag_summary.tsv")
        proc = subprocess.run([sys.executable, SCRIPT, "-b", self.bam, "-g", gtf, "-o", out["bam"],
                               "-s", out["tsv"], "-t", "1", *extra_args], capture_output=True, text=True)
        return proc, out

    def read_tsv(self, path):
        with open(path) as fh:
            return [line.rstrip("\n").split("\t") for line in fh]

    def tags(self, bam):
        with pysam.AlignmentFile(bam, check_sq=False) as fh:
            return {r.query_name: (r.get_tag("GX"), r.get_tag("GS")) for r in fh.fetch(until_eof=True)}

    def test_matching_names(self):
        proc, out = self.run_tool("match", GENES)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertNotIn("WARNING", proc.stderr)
        self.assertIn("contig overlap: 1 of 1 gtf seqnames", proc.stderr)

        tags = self.tags(out["bam"])
        self.assertEqual({name: tags[name][1] for name in tags}, {r[0]: r[4] for r in READS})
        self.assertEqual(tags["in_g1"][0], "G1")
        self.assertEqual(tags["in_g2_g3"][0], "G2")
        self.assertEqual(tags["intergenic"][0], "__no_gene")

        # 2 of the 4 alignments carry no gene
        self.assertIn("GS=none fraction: 0.5000", proc.stderr)

        with open(out["tsv"]) as fh:
            rows = [line.rstrip("\n").split("\t") for line in fh]
        self.assertEqual(rows[0], ["contig", "gs", "non_primary", "count"])
        self.assertIn(["chr2", "none", "0", "1"], rows[1:])

    def test_no_tag_summary_by_default(self):
        proc, out = self.run_tool("nosummary", GENES)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertFalse(os.path.exists(out["tag_summary"]))

    def test_tag_summary(self):
        tag_summary = os.path.join(self.tmp.name, "S1.cdna.genome.gene_tag_summary.tsv")
        proc, _ = self.run_tool("summary", GENES, ["-q", tag_summary, "-n", "S1"])
        self.assertEqual(proc.returncode, 0, proc.stderr)

        rows = self.read_tsv(tag_summary)
        self.assertEqual(rows[0], ["sample", "reads", "gs_none_reads", "gs_none_fraction"])
        # 2 of the 4 alignments carry no gene, the same share the log reports
        self.assertEqual(rows[1:], [["S1", "4", "2", "0.500000"]])
        self.assertIn("GS=none fraction: 0.5000", proc.stderr)

    def test_tag_summary_sample_defaults_to_file_name(self):
        tag_summary = os.path.join(self.tmp.name, "S2.a.cdna.genome.gene_tag_summary.tsv")
        proc, _ = self.run_tool("default_name", GENES, ["--tag_summary", tag_summary])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.read_tsv(tag_summary)[1][0], "S2.a.cdna.genome")

    def test_tag_summary_without_genes_on_any_read(self):
        # every read is intergenic for this gtf: the whole sample is GS=none
        far = [("chr1", 90001, 91000, "GF", "GeneFar")]
        tag_summary = os.path.join(self.tmp.name, "far.gene_tag_summary.tsv")
        proc, _ = self.run_tool("far", far, ["-q", tag_summary, "-n", "far"])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.read_tsv(tag_summary)[1], ["far", "4", "4", "1.000000"])

    def test_chr_prefix_mismatch_fails(self):
        # Ensembl-style gtf ("1") against a UCSC-style bam ("chr1"): no shared contig
        tag_summary = os.path.join(self.tmp.name, "ensembl.gene_tag_summary.tsv")
        proc, out = self.run_tool("ensembl", renamed(GENES, lambda c: c[3:]), ["-q", tag_summary])
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("none of the 1 gtf seqnames", proc.stderr)
        self.assertIn("'chr' prefix", proc.stderr)
        # it fails before writing anything
        self.assertFalse(os.path.exists(out["bam"]))
        self.assertFalse(os.path.exists(out["tsv"]))
        self.assertFalse(os.path.exists(out["tag_summary"]))

    def test_small_overlap_warns(self):
        # the gtf annotates chr1 plus three contigs the bam does not have: 1 of 4 shared
        extra = [(contig, 100, 200, f"X{contig}", f"X{contig}") for contig in ("chr3", "chr4", "chr5")]
        proc, out = self.run_tool("partial", GENES + extra)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("WARNING: only 25.0% of the gtf seqnames", proc.stderr)
        self.assertIn("chr3", proc.stderr)
        self.assertTrue(os.path.exists(out["bam"]))

    def test_half_overlap_does_not_warn(self):
        # the threshold is a strict "below half"
        proc, _ = self.run_tool("half", GENES + [("chr9", 100, 200, "X9", "X9")])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("(50.0%)", proc.stderr)
        self.assertNotIn("WARNING", proc.stderr)


if __name__ == "__main__":
    unittest.main()

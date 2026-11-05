#!/usr/bin/env python3

""" Unit test for bin/pre_extract_barcodes.py (BLAZE path) on toy reads.

    Not run by nf-test (which only collects *.nf.test); .github/workflows/bin-tests.yml
    runs it on pull requests. Locally, run it in the module's container, which has
    biopython but no pytest:

        apptainer exec <mulled-v2-497e12343495a3e0f3b3459542cc8ad37575d9fa image> \
            python -m unittest tests/bin/test_pre_extract_barcodes.py
"""

import os
import subprocess
import sys
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin", "pre_extract_barcodes.py")

READ1 = "CTACACGACGCTCTTCCGATCT"
TSO = "TTTCTTATATGGG"
BC = "AACCGGTTACGTCCAA"     # not its own reverse complement
CDNA = "GATTACAGATTACAGATTACAGATTACA"
RC = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(RC)[::-1]


def write_inputs(d, reads, csv_ids=None):
    """reads: (read id, sequence, barcode); csv_ids overrides the ids in the BLAZE csv."""
    fq = os.path.join(d, "in.fastq")
    csv = os.path.join(d, "in.csv")
    with open(fq, "w") as fh:
        for read_id, seq, _ in reads:
            fh.write(f"@{read_id} extra=1\n{seq}\n+\n{'5' * len(seq)}\n")
    ids = csv_ids if csv_ids is not None else [r[0] for r in reads]
    with open(csv, "w") as fh:
        for read_id, (_, _, bc) in zip(ids, reads):
            # BLAZE: read_id,putative_bc,putative_bc_min_q,putative_umi,polyT_end,pre_bc_flanking,post_umi_flanking
            fh.write(f"{read_id},{bc},20,,,,\n")
    return fq, csv


class PreExtractBarcodesTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def run_tool(self, fq, csv, fmt, check=True):
        prefix = os.path.join(self.d, "out")
        res = subprocess.run([sys.executable, SCRIPT, "-i", fq, "-b", csv, "-o", prefix, "-f", fmt, "-t", "1"],
                             capture_output=True, text=True)
        if check:
            self.assertEqual(res.returncode, 0, res.stderr)
        return prefix, res

    def outputs(self, prefix):
        with open(prefix + ".putative_bc_umi.tsv") as fh:
            tsv = [line.rstrip("\n").split("\t") for line in fh][1:]
        with open(prefix + ".fastq") as fh:
            lines = fh.read().splitlines()
        fastq = {lines[i][1:]: lines[i + 1] for i in range(0, len(lines), 4)}
        return {row[0]: row for row in tsv}, fastq

    def test_5v3_umi_and_tso(self):
        umi = "ACGTACGTAC"
        reads = [("fwd", READ1 + BC + umi + TSO + CDNA, BC),
                 ("rev", revcomp(READ1 + BC + umi + TSO + CDNA), BC),
                 ("nobc", READ1 + CDNA, "")]
        prefix, _ = self.run_tool(*write_inputs(self.d, reads), "10X_5v3")
        tsv, fastq = self.outputs(prefix)
        self.assertEqual(set(tsv), {"fwd", "rev"})
        for read_id in ("fwd", "rev"):
            self.assertEqual(tsv[read_id][1], BC)
            # a 10 bp UMI, then the whole 13 bp TSO is stripped
            self.assertEqual(tsv[read_id][3], umi)
            self.assertEqual(fastq[read_id], CDNA)

    def test_5v2_matches_5v3(self):
        umi = "ACGTACGTAC"
        reads = [("fwd", READ1 + BC + umi + TSO + CDNA, BC)]
        prefix, _ = self.run_tool(*write_inputs(self.d, reads), "10X_5v2")
        tsv, fastq = self.outputs(prefix)
        self.assertEqual(tsv["fwd"][3], umi)
        self.assertEqual(fastq["fwd"], CDNA)

    def test_3v3_umi_and_polyt(self):
        umi = "ACGTACGTACGT"
        reads = [("fwd", READ1 + BC + umi + "T" * 10 + CDNA, BC)]
        prefix, _ = self.run_tool(*write_inputs(self.d, reads), "10X_3v3")
        tsv, fastq = self.outputs(prefix)
        self.assertEqual(tsv["fwd"][3], umi)
        self.assertEqual(fastq["fwd"], CDNA)

    def test_read_id_mismatch_fails(self):
        reads = [("a", READ1 + BC + CDNA, BC), ("b", READ1 + BC + CDNA, BC)]
        _, res = self.run_tool(*write_inputs(self.d, reads, csv_ids=["a", "c"]), "10X_3v3", check=False)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("read 2 is b", res.stderr)

    def test_read_count_mismatch_fails(self):
        reads = [("a", READ1 + BC + CDNA, BC), ("b", READ1 + BC + CDNA, BC)]
        _, res = self.run_tool(*write_inputs(self.d, reads, csv_ids=["a"]), "10X_3v3", check=False)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("different number of reads", res.stderr)

    def test_unknown_format_rejected(self):
        # an unknown format used to drop every read silently
        reads = [("a", READ1 + BC + CDNA, BC)]
        _, res = self.run_tool(*write_inputs(self.d, reads), "cellranger", check=False)
        self.assertNotEqual(res.returncode, 0)


if __name__ == "__main__":
    unittest.main()

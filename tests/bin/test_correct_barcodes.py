#!/usr/bin/env python3

""" Unit test for bin/correct_barcodes.py (BLAZE path) on a toy putative barcode table.

    Not run by nf-test (which only collects *.nf.test); .github/workflows/bin-tests.yml
    runs it on pull requests. Locally, run it in the module's container, which has pygtrie
    but no pytest:

        apptainer exec <mulled-v2-bb96c7354781ab52d8e69ccff89587598dc87fea image> \
            python -m unittest tests/bin/test_correct_barcodes.py
"""

import os
import subprocess
import sys
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin", "correct_barcodes.py")

HEADER = "read_id\tbc\tbc_qual\tumi\tumi_qual"
# whitelisted barcodes and their counts; ACGTACGA and ACGTACGC are one substitution apart
# from ACGTACGG, and ACGTACGA is three times as abundant
WHITELIST = {"ACGTACGA": 3, "ACGTACGC": 1, "TTTTGGGG": 5}
QUAL = "+" * 8      # Q10 everywhere
# read id, barcode, expected correction with the default --min_post_prob (None: dropped)
READS = [
    ("exact", "TTTTGGGG", "TTTTGGGG"),
    ("one_mismatch", "TTTTGGGC", "TTTTGGGG"),
    # posterior 4/6 for ACGTACGA vs 2/6 for ACGTACGC: below 0.975
    ("ambiguous", "ACGTACGG", None),
    ("far", "CCCCCCCC", None),
    ("no_bc", "", None),
]


def write_inputs(d, header=True):
    paths = {k: os.path.join(d, k) for k in ("in.tsv", "whitelist.txt", "bc_count.txt")}
    with open(paths["in.tsv"], "w") as fh:
        if header:
            fh.write(HEADER + "\n")
        for read_id, bc, _ in READS:
            fh.write("\t".join([read_id, bc, QUAL if bc else "", "AAAAAAAAAAAA", "I" * 12]) + "\n")
    with open(paths["whitelist.txt"], "w") as fh:
        fh.write("".join(bc + "\n" for bc in WHITELIST))
    with open(paths["bc_count.txt"], "w") as fh:
        fh.write("".join(f"{bc},{n}\n" for bc, n in WHITELIST.items()))
    return paths


class CorrectBarcodesTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def run_tool(self, paths, *extra):
        out = os.path.join(self.d, "out.tsv")
        res = subprocess.run([sys.executable, SCRIPT, "-i", paths["in.tsv"], "-o", out,
                              "-w", paths["whitelist.txt"], "-b", paths["bc_count.txt"], *extra],
                             check=True, capture_output=True, text=True)
        with open(out) as fh:
            lines = [line.rstrip("\n").split("\t") for line in fh]
        return lines, res.stdout

    def test_corrections(self):
        lines, _ = self.run_tool(write_inputs(self.d), "--skip_header")
        got = {line[0]: line[-1] for line in lines}
        self.assertEqual(got, {r: c for r, _, c in READS if c})

    def test_decimal_min_post_prob(self):
        # was type=int, so any decimal on the command line failed argparse
        lines, _ = self.run_tool(write_inputs(self.d), "--skip_header", "--min_post_prob", "0.6")
        got = {line[0]: line[-1] for line in lines}
        self.assertEqual(got["ambiguous"], "ACGTACGA")

    def test_header_written_unless_skipped(self):
        lines, _ = self.run_tool(write_inputs(self.d))
        self.assertEqual(lines[0], HEADER.split("\t") + ["corrected_bc"])
        self.assertEqual(len(lines), 3)
        lines, _ = self.run_tool(write_inputs(self.d), "--skip_header")
        self.assertEqual(len(lines), 2)

    def test_input_without_header_keeps_first_read(self):
        # the first line used to be dropped as a header whatever it was
        lines, _ = self.run_tool(write_inputs(self.d, header=False), "--skip_header")
        self.assertEqual([line[0] for line in lines], ["exact", "one_mismatch"])

    def test_one_summary_line(self):
        # no line per candidate barcode
        _, stdout = self.run_tool(write_inputs(self.d), "--skip_header")
        self.assertEqual(stdout.splitlines(), [
            "INFO: 5 reads: 1 on the whitelist, 1 corrected, 2 not corrected (dropped), "
            "1 without a barcode (dropped)"])


if __name__ == "__main__":
    unittest.main()

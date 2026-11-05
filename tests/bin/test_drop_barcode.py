#!/usr/bin/env python3

""" Unit test for bin/drop_barcode.py on small MEX and tsv fixtures with a '-' barcode.

    Not run by nf-test (which only collects *.nf.test); .github/workflows/bin-tests.yml
    runs it on pull requests. It needs only the python standard library; locally, run it
    in the module's container (modules/local/drop_barcode):

        apptainer exec <python:3.10.4 image> python -m unittest -v tests/bin/test_drop_barcode.py

    The fixtures follow IsoQuant 3.13 `--counts_format mtx` output (features on the rows,
    barcodes on the columns, flexiplex's '-' among the barcodes), its transposed layout,
    gzipped copies of both, and the `linear` and wide tsv tables.
"""

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin", "drop_barcode.py")

FEATURES = ["ENSG01", "ENSG02", "ENSG03"]
# '-' second, so the barcodes after it have to be renumbered
BARCODES = ["AAAA", "-", "CCCC", "GGGG"]
# (feature index, barcode index, value), 1-based, as IsoQuant writes them
ENTRIES = [(1, 1, "1.0"), (1, 2, "5.0"), (2, 2, "2.0"), (2, 3, "1.0"),
           (3, 4, "3.0"), (3, 2, "0.5"), (2, 1, "4.0")]
BANNER = "%%MatrixMarket matrix coordinate real general\n"


def mtx_text(entries, nrows, ncols, transpose=False):
    lines = [BANNER, "% written by a test\n"]
    if transpose:
        nrows, ncols = ncols, nrows
        entries = [(j, i, v) for i, j, v in entries]
    lines.append(f"{nrows} {ncols} {len(entries)}\n")
    lines += [f"{i} {j} {v}\n" for i, j, v in entries]
    return "".join(lines)


def read_text(path):
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode()


def read_mtx(path, transpose=False):
    """Return (banner and comments, dims, {(feature_name_idx, barcode_idx): value})"""
    lines = read_text(path).splitlines()
    header = [l for l in lines if l.startswith("%")]
    body = [l for l in lines if not l.startswith("%")]
    dims = tuple(int(x) for x in body[0].split())
    cells = {}
    for line in body[1:]:
        i, j, v = line.split()
        key = (int(j), int(i)) if transpose else (int(i), int(j))
        cells[key] = float(v)
    return header, dims, cells


def dense_by_name(mtx, barcodes, features, transpose=False):
    """{(feature, barcode): value}, so outputs can be compared independently of indices"""
    _, _, cells = read_mtx(mtx, transpose)
    bcs = read_text(barcodes).splitlines()
    feats = [l.split("\t")[0] for l in read_text(features).splitlines()]
    return {(feats[i - 1], bcs[j - 1]): v for (i, j), v in cells.items()}


class DropBarcodeTest(unittest.TestCase):

    def setUp(self):
        self.workdir = tempfile.mkdtemp(prefix="drop_barcode_", dir=os.environ.get("TMPDIR"))
        self.addCleanup(shutil.rmtree, self.workdir)
        self.indir = os.path.join(self.workdir, "in")
        self.outdir = os.path.join(self.workdir, "out")
        os.makedirs(self.indir)

    def write(self, name, text, gz=False):
        path = os.path.join(self.indir, name)
        data = text.encode()
        with open(path, "wb") as fh:
            fh.write(gzip.compress(data, mtime=0) if gz else data)
        return path

    def mex(self, barcodes=BARCODES, entries=ENTRIES, transpose=False, gz=False):
        ext = ".gz" if gz else ""
        nrows, ncols = len(FEATURES), len(barcodes)
        return (
            self.write(f"s.gene_grouped_tag_CB_counts.matrix.mtx{ext}", mtx_text(entries, nrows, ncols, transpose), gz),
            self.write(f"s.gene_grouped_tag_CB_counts.barcodes.tsv{ext}", "".join(b + "\n" for b in barcodes), gz),
            self.write(f"s.gene_grouped_tag_CB_counts.features.tsv{ext}", "".join(f"{f}\t{f}\n" for f in FEATURES), gz),
        )

    def run_script(self, *args, check=True):
        proc = subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)
        if check:
            self.assertEqual(proc.returncode, 0, proc.stderr)
        return proc

    def run_mex(self, files, *extra, check=True):
        mtx, bcs, feats = files
        return self.run_script("mex", "-m", mtx, "-b", bcs, "-f", feats, "-o", self.outdir, *extra, check=check)

    def out(self, path):
        return os.path.join(self.outdir, os.path.basename(path))

    def assert_dash_dropped(self, files, transpose=False):
        mtx, bcs, feats = files
        before = dense_by_name(mtx, bcs, feats, transpose)
        after = dense_by_name(self.out(mtx), self.out(bcs), self.out(feats), transpose)
        # every count on a real barcode is kept, and nothing else
        self.assertEqual(after, {k: v for k, v in before.items() if k[1] != "-"})
        self.assertEqual(read_text(self.out(bcs)).splitlines(), ["AAAA", "CCCC", "GGGG"])
        with open(feats, "rb") as a, open(self.out(feats), "rb") as b:
            self.assertEqual(a.read(), b.read())
        header, dims, cells = read_mtx(self.out(mtx), transpose)
        self.assertEqual(header, [BANNER.strip(), "% written by a test"])
        self.assertEqual(dims, (3, 3, 4))
        # the remaining barcodes are renumbered 1..3 without gaps
        self.assertEqual(sorted({j for _, j in cells}), [1, 2, 3])
        self.assertEqual(cells[(3, 3)], 3.0)

    def test_barcodes_on_columns(self):
        files = self.mex()
        proc = self.run_mex(files)
        self.assert_dash_dropped(files)
        self.assertIn("col 2 of 4", proc.stderr)
        self.assertIn("3 of 7 entries, 7.5 counts", proc.stderr)

    def test_barcodes_on_rows(self):
        files = self.mex(transpose=True)
        proc = self.run_mex(files)
        self.assert_dash_dropped(files, transpose=True)
        self.assertIn("row 2 of 4", proc.stderr)

    def test_gzipped(self):
        files = self.mex(gz=True)
        self.run_mex(files)
        for path in files:
            with open(self.out(path), "rb") as fh:
                self.assertEqual(fh.read(2), b"\x1f\x8b", path)
        self.assert_dash_dropped(files)

    def test_gzipped_without_gz_suffix(self):
        files = self.mex(gz=True)
        renamed = []
        for path in files:
            shutil.move(path, path[:-3])
            renamed.append(path[:-3])
        self.run_mex(renamed)
        self.assert_dash_dropped(renamed)

    def test_dash_first_and_last(self):
        barcodes = ["-", "AAAA", "CCCC", "GGGG", "-"]
        entries = [(1, 1, "2"), (1, 2, "1"), (2, 5, "7"), (3, 4, "3"), (2, 3, "1")]
        files = self.mex(barcodes=barcodes, entries=entries)
        self.run_mex(files)
        mtx, bcs, feats = files
        _, dims, cells = read_mtx(self.out(mtx))
        self.assertEqual(dims, (3, 3, 3))
        self.assertEqual(cells, {(1, 1): 1.0, (3, 3): 3.0, (2, 2): 1.0})
        self.assertEqual(read_text(self.out(bcs)).splitlines(), ["AAAA", "CCCC", "GGGG"])

    def test_no_dash_is_copied_unchanged(self):
        for gz in (False, True):
            with self.subTest(gz=gz):
                barcodes = ["AAAA", "CCCC", "GGGG", "TTTT"]
                files = self.mex(barcodes=barcodes, gz=gz)
                proc = self.run_mex(files)
                self.assertIn("copied unchanged", proc.stderr)
                for path in files:
                    with open(path, "rb") as a, open(self.out(path), "rb") as b:
                        self.assertEqual(a.read(), b.read(), path)

    def test_square_matrix(self):
        # 3 features x 3 barcodes: the dimensions cannot tell which axis holds the barcodes
        barcodes = ["AAAA", "-", "CCCC"]
        entries = [(1, 2, "1"), (2, 1, "1"), (3, 3, "1")]
        files = self.mex(barcodes=barcodes, entries=entries)
        proc = self.run_mex(files)
        self.assertIn("either dimension", proc.stderr)
        _, dims, cells = read_mtx(self.out(files[0]))
        self.assertEqual((dims, cells), ((3, 2, 2), {(2, 1): 1.0, (3, 2): 1.0}))

        shutil.rmtree(self.outdir)
        self.run_mex(files, "--barcode_axis", "rows")
        _, dims, cells = read_mtx(self.out(files[0]))
        self.assertEqual((dims, cells), ((2, 3, 2), {(1, 2): 1.0, (2, 3): 1.0}))

    def test_without_features(self):
        mtx, bcs, _ = self.mex()
        self.run_script("mex", "-m", mtx, "-b", bcs, "-o", self.outdir)
        _, dims, _ = read_mtx(self.out(mtx))
        self.assertEqual(dims, (3, 3, 4))

    def test_other_barcode(self):
        files = self.mex()
        self.run_mex(files, "--barcode", "CCCC")
        self.assertEqual(read_text(self.out(files[1])).splitlines(), ["AAAA", "-", "GGGG"])

    def test_dimension_mismatch_fails(self):
        files = self.mex()
        self.write(os.path.basename(files[1]), "AAAA\n-\nCCCC\n")
        proc = self.run_mex(files, check=False)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("does not match 3 barcodes", proc.stderr)

    def test_entry_count_mismatch_fails(self):
        files = self.mex()
        text = read_text(files[0]).replace("3 4 7\n", "3 4 8\n")
        self.write(os.path.basename(files[0]), text)
        proc = self.run_mex(files, check=False)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("says 8 entries but the file has 7", proc.stderr)

    def test_refuses_to_overwrite_input(self):
        mtx, bcs, feats = self.mex()
        proc = self.run_script("mex", "-m", mtx, "-b", bcs, "-f", feats, "-o", self.indir, check=False)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("would overwrite its input", proc.stderr)

    #
    # tsv
    #

    def test_linear_tsv(self):
        text = ("feature_id\tgroup_id\tcount\n"
                "ENSG01\tAAAA\t1.00\n"
                "ENSG01\t-\t5.00\n"
                "ENSG02\tCCCC\t1.00\n"
                "ENSG03\t-\t0.50\n")
        path = self.write("s.gene_grouped_tag_CB_counts.linear.tsv", text)
        proc = self.run_script("tsv", "-i", path, "-o", self.outdir)
        self.assertIn("2 row(s)", proc.stderr)
        self.assertEqual(read_text(self.out(path)),
                         "feature_id\tgroup_id\tcount\nENSG01\tAAAA\t1.00\nENSG02\tCCCC\t1.00\n")

    def test_wide_tsv_barcode_columns(self):
        text = ("#feature_id\tAAAA\t-\tCCCC\n"
                "ENSG01\t1\t5\t0\n"
                "ENSG02\t0\t2\t1\n")
        path = self.write("s.gene_grouped_tag_CB_counts.tsv", text, gz=True)
        self.run_script("tsv", "-i", path, "-o", self.outdir)
        with open(self.out(path), "rb") as fh:
            self.assertEqual(fh.read(2), b"\x1f\x8b")
        self.assertEqual(read_text(self.out(path)), "#feature_id\tAAAA\tCCCC\nENSG01\t1\t0\nENSG02\t0\t1\n")

    def test_wide_tsv_barcode_rows(self):
        text = ("barcode_id\tENSG01\tENSG02\n"
                "AAAA\t1\t0\n"
                "-\t5\t2\n"
                "CCCC\t0\t1\n")
        path = self.write("s.tsv", text)
        self.run_script("tsv", "-i", path, "-o", self.outdir)
        self.assertEqual(read_text(self.out(path)), "barcode_id\tENSG01\tENSG02\nAAAA\t1\t0\nCCCC\t0\t1\n")

    def test_tsv_without_dash_is_copied_unchanged(self):
        for name, text in (("linear.tsv", "feature_id\tgroup_id\tcount\nENSG01\tAAAA\t1.00\n"),
                           ("wide.tsv", "#feature_id\tAAAA\tCCCC\r\nENSG01\t1\t0\r\n")):
            with self.subTest(name=name):
                path = self.write(name, text)
                proc = self.run_script("tsv", "-i", path, "-o", self.outdir)
                self.assertIn("copied unchanged", proc.stderr)
                with open(path, "rb") as a, open(self.out(path), "rb") as b:
                    self.assertEqual(a.read(), b.read())


if __name__ == "__main__":
    unittest.main()

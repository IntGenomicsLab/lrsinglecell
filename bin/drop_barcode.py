#!/usr/bin/env python3

""" Remove one barcode from a feature-barcode count matrix.

    Written for flexiplex's '-' placeholder: with `-a true` flexiplex leaves `CB` as '-'
    on every read that matched no known barcode, so a count matrix grouped on `CB` gets
    one '-' pseudo-cell holding all of them. This drops that barcode and keeps every
    other count as it was.

    Two input kinds are handled, plain or gzipped (detected from the content, not the
    file name); outputs are written under --outdir with the input file names and the
    input compression:

      mex  A Matrix Market coordinate matrix with its barcodes and (optionally) features
           list. Which matrix dimension holds the barcodes is worked out from the
           dimensions in the matrix header against the length of the barcodes (and
           features) list, not assumed. The barcode's row or column is removed, the
           remaining indices are renumbered to stay contiguous, the header's dimension
           and entry counts are updated and the barcode is removed from the list. The
           features list is copied unchanged.

      tsv  An IsoQuant-style counts table. A `linear` table (a header naming a
           `group_id` or `barcode` column) loses the rows for the barcode. A wide table
           loses the column named after the barcode or, if the barcodes are the rows,
           the row whose first field is the barcode.

    A file that does not contain the barcode is copied byte for byte.
"""

import argparse
import gzip
import io
import os
import shutil
import sys

GZIP_MAGIC = b"\x1f\x8b"
LINEAR_BARCODE_COLUMNS = ("group_id", "barcode", "barcodes", "cell", "cell_id")


def parse_args(argv=None):
    """Parse the commandline arguments"""

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-c", "--barcode", default="-",
                        help="The barcode to remove (default: '-', flexiplex's no-match placeholder)")

    arg_parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    sub = arg_parser.add_subparsers(dest="kind", required=True)

    mex = sub.add_parser("mex", parents=[common], help="A Matrix Market matrix with its barcodes and features lists")
    mex.add_argument("-m", "--matrix", required=True, help="The matrix.mtx(.gz)")
    mex.add_argument("-b", "--barcodes", required=True, help="The barcodes.tsv(.gz), one barcode per line")
    mex.add_argument("-f", "--features", default=None, help="The features.tsv(.gz); used to check dimensions and copied")
    mex.add_argument("-a", "--barcode_axis", choices=("auto", "rows", "cols"), default="auto",
                     help="The matrix dimension that holds the barcodes (default: worked out from the files)")
    mex.add_argument("-o", "--outdir", required=True, help="Output directory")

    tsv = sub.add_parser("tsv", parents=[common], help="An IsoQuant linear or wide counts table")
    tsv.add_argument("-i", "--input", required=True, help="The counts table")
    tsv.add_argument("-o", "--outdir", required=True, help="Output directory")

    return arg_parser.parse_args(argv)


def fail(message):
    sys.exit(f"drop_barcode.py: error: {message}")


def log(message):
    print(f"drop_barcode.py: {message}", file=sys.stderr)


def is_gzipped(path):
    with open(path, "rb") as fh:
        return fh.read(2) == GZIP_MAGIC


def open_read(path):
    """Text reader that keeps line endings as they are"""
    if is_gzipped(path):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


class open_write:
    """Text writer with the compression of `like`; gzip output carries no name or mtime"""

    def __init__(self, path, like):
        self.path = path
        self.gz = is_gzipped(like)

    def __enter__(self):
        self.raw = open(self.path, "wb")
        stream = gzip.GzipFile(filename="", mode="wb", fileobj=self.raw, mtime=0) if self.gz else self.raw
        self.fh = io.TextIOWrapper(stream, encoding="utf-8", newline="")
        return self.fh

    def __exit__(self, *exc):
        self.fh.close()
        if not self.raw.closed:
            self.raw.close()
        return False


def output_path(outdir, path):
    out = os.path.join(outdir, os.path.basename(path))
    if os.path.exists(out) and os.path.samefile(out, path):
        fail(f"output {out} would overwrite its input; choose another --outdir")
    return out


def split_ending(line):
    body = line.rstrip("\r\n")
    return body, line[len(body):]


def copy_unchanged(path, outdir):
    out = output_path(outdir, path)
    shutil.copyfile(path, out)
    return out


#
# MEX
#

def first_fields(path):
    """The first tab-separated field of every line, ignoring blank lines at the end"""
    with open_read(path) as fh:
        fields = [split_ending(line)[0].split("\t")[0] for line in fh]
    while fields and not fields[-1].strip():
        fields.pop()
    return fields


def read_mtx_header(fh):
    """Return the header lines (banner, comments, size line) and the parsed size line"""
    banner = fh.readline()
    tokens = banner.lower().split()
    if len(tokens) < 5 or tokens[0] != "%%matrixmarket" or tokens[1] != "matrix":
        fail("matrix is not a Matrix Market file (no '%%MatrixMarket matrix' banner)")
    if tokens[2] != "coordinate":
        fail(f"only coordinate Matrix Market files are supported, not '{tokens[2]}'")
    if tokens[4] != "general":
        fail(f"only 'general' matrices can lose a barcode, not '{tokens[4]}'")
    header = [banner]
    for line in fh:
        header.append(line)
        if not line.startswith("%") and line.strip():
            try:
                nrows, ncols, nnz = (int(x) for x in line.split())
            except ValueError:
                fail(f"malformed Matrix Market size line: {line.strip()!r}")
            return header, (nrows, ncols, nnz)
    fail("matrix has no size line")


def barcode_axis(nrows, ncols, n_barcodes, n_features, requested):
    """Which matrix dimension ('rows' or 'cols') the barcodes list indexes"""
    dims = {"rows": nrows, "cols": ncols}
    other = {"rows": "cols", "cols": "rows"}
    candidates = [requested] if requested != "auto" else ["cols", "rows"]
    fits = [axis for axis in candidates
            if dims[axis] == n_barcodes and (n_features is None or dims[other[axis]] == n_features)]
    if not fits:
        described = f"{n_barcodes} barcodes" + ("" if n_features is None else f" and {n_features} features")
        where = "" if requested == "auto" else f" with the barcodes on the {requested}"
        fail(f"a {nrows} x {ncols} matrix does not match {described}{where}")
    if len(fits) > 1:
        # square, and the features list (if any) cannot tell the two apart
        log(f"a {nrows} x {ncols} matrix fits the barcodes on either dimension; "
            "taking the columns (the MEX convention). Use --barcode_axis to choose")
    return fits[0]


def drop_from_mex(matrix, barcodes, features, outdir, barcode, requested_axis="auto"):
    """Write the matrix, barcodes and features without `barcode`; return a summary dict"""
    bc_list = first_fields(barcodes)
    drop = [i + 1 for i, bc in enumerate(bc_list) if bc == barcode]

    out_matrix = output_path(outdir, matrix)
    out_barcodes = output_path(outdir, barcodes)
    out_features = output_path(outdir, features) if features else None

    if not drop:
        for path in (matrix, barcodes, features):
            if path:
                copy_unchanged(path, outdir)
        log(f"{os.path.basename(barcodes)}: no '{barcode}' barcode, copied unchanged")
        return {"dropped_barcodes": 0, "dropped_entries": 0, "dropped_total": 0.0}

    n_features = len(first_fields(features)) if features else None
    with open_read(matrix) as fh:
        _, (nrows, ncols, nnz) = read_mtx_header(fh)
    axis = barcode_axis(nrows, ncols, len(bc_list), n_features, requested_axis)
    col = 1 if axis == "cols" else 0

    # old 1-based index -> new 1-based index, 0 for a dropped barcode
    remap = [0] * (len(bc_list) + 1)
    dropped = set(drop)
    new = 0
    for old in range(1, len(bc_list) + 1):
        if old not in dropped:
            new += 1
            remap[old] = new

    # Two passes over the entries, so the size line can be written before them
    dropped_entries = 0
    dropped_total = 0.0
    seen = 0
    with open_read(matrix) as fh:
        read_mtx_header(fh)
        for line in fh:
            fields = line.split()
            if not fields or fields[0].startswith("%"):
                continue
            seen += 1
            if remap[int(fields[col])] == 0:
                dropped_entries += 1
                dropped_total += float(fields[2]) if len(fields) > 2 else 1.0
    if seen != nnz:
        fail(f"matrix header says {nnz} entries but the file has {seen}")

    with open_read(matrix) as fh, open_write(out_matrix, matrix) as out:
        header, _ = read_mtx_header(fh)
        size_end = split_ending(header[-1])[1] or "\n"
        dims = [nrows, ncols]
        dims[col] -= len(drop)
        out.writelines(header[:-1])
        out.write(f"{dims[0]} {dims[1]} {nnz - dropped_entries}{size_end}")
        for line in fh:
            body, end = split_ending(line)
            fields = body.split()
            if not fields or fields[0].startswith("%"):
                out.write(line)
                continue
            idx = remap[int(fields[col])]
            if idx == 0:
                continue
            fields[col] = str(idx)
            out.write(" ".join(fields) + end)

    with open_read(barcodes) as fh, open_write(out_barcodes, barcodes) as out:
        for line in fh:
            if split_ending(line)[0].split("\t")[0] != barcode:
                out.write(line)

    if features:
        copy_unchanged(features, outdir)

    log(f"{os.path.basename(matrix)}: dropped barcode '{barcode}' ({axis[:-1]} "
        f"{', '.join(map(str, drop))} of {len(bc_list)}): {dropped_entries} of {nnz} entries, "
        f"{dropped_total:g} counts; {len(bc_list) - len(drop)} barcodes kept")
    return {"dropped_barcodes": len(drop), "dropped_entries": dropped_entries, "dropped_total": dropped_total,
            "axis": axis, "outputs": [p for p in (out_matrix, out_barcodes, out_features) if p]}


#
# TSV
#

def drop_from_tsv(table, outdir, barcode):
    """Write the table without `barcode`; return a summary dict"""
    out_table = output_path(outdir, table)

    with open_read(table) as fh:
        header = fh.readline()
        names = split_ending(header)[0].split("\t")
        keys = [n.lstrip("#").strip().lower() for n in names]
        linear = next((keys.index(k) for k in LINEAR_BARCODE_COLUMNS if k in keys), None)
        if linear is not None:
            mode, drop_cols = "linear", []
            present = any(split_ending(line)[0].split("\t")[linear:linear + 1] == [barcode] for line in fh)
        else:
            drop_cols = [i for i, n in enumerate(names) if i > 0 and n == barcode]
            in_rows = any(split_ending(line)[0].split("\t")[0] == barcode for line in fh)
            if drop_cols and in_rows:
                fail(f"{table}: '{barcode}' is both a column name and a row name; cannot tell which holds the barcodes")
            mode = "columns" if drop_cols else "rows"
            present = bool(drop_cols) or in_rows

    if not present:
        copy_unchanged(table, outdir)
        log(f"{os.path.basename(table)}: no '{barcode}' barcode, copied unchanged")
        return {"dropped_barcodes": 0, "mode": mode}

    dropped_rows = 0
    key = linear if mode == "linear" else 0
    with open_read(table) as fh, open_write(out_table, table) as out:
        for n, line in enumerate(fh):
            body, end = split_ending(line)
            fields = body.split("\t")
            if mode == "columns":
                out.write("\t".join(f for i, f in enumerate(fields) if i not in drop_cols) + end)
            elif n > 0 and fields[key:key + 1] == [barcode]:
                dropped_rows += 1
            else:
                out.write(line)

    what = f"{len(drop_cols)} column(s)" if mode == "columns" else f"{dropped_rows} row(s)"
    log(f"{os.path.basename(table)}: dropped barcode '{barcode}' ({mode} table): {what}")
    return {"dropped_barcodes": 1, "mode": mode, "dropped_rows": dropped_rows, "outputs": [out_table]}


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    if args.kind == "mex":
        drop_from_mex(args.matrix, args.barcodes, args.features, args.outdir, args.barcode, args.barcode_axis)
    else:
        drop_from_tsv(args.input, args.outdir, args.barcode)


if __name__ == "__main__":
    main()

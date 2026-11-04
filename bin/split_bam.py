#!/usr/bin/env python3

# Written by Austyn Trull and released under the MIT license.

import argparse
import pysam

def parse_args():
    """ Parse commandline arguments"""

    parser = argparse.ArgumentParser(description="Split a BAM file into multiple BAM files based on read groups.")

    parser.add_argument(
        "-i",
        "--input",
        default=None,
        type=str,
        required=True,
        help="Input BAM file"
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        type=str,
        required=True,
        help="Output BAM file"
    )

    parser.add_argument(
        "-c",
        "--contig_file",
        default=None,
        type=str,
        required=True,
        help="File containing contig names to include in the output BAM file."
    )

    parser.add_argument(
        "--counts",
        default=None,
        type=str,
        required=False,
        help="Optional file to write '<input reads>\\t<reads written>' to, so the caller can "
             "check that the chunks of one BAM are not all empty."
    )

    return parser.parse_args()

def split_bam_by_contig(input_bam, output_bam, contigs_file):
    """ Split a BAM file into multiple BAM files based on contig names.

    Args:
        input_bam (str): Path to the input BAM file.
        output_bam (str): Path to the output BAM file.
        contigs_file (str): Path to the file containing contig names to include in the output BAM file.

    Returns:
        tuple: (number of reads in the input BAM, number of reads written)
    """

    # Read the contig names from the contigs file
    contigs_to_include = set()
    with open(contigs_file, 'r') as f:
        contigs_to_include = set(line.strip() for line in f)

    # Open the input BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")

    # Create a new BAM file for output
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    n_input = 0
    n_written = 0

    # Iterate through each read in the input BAM file
    for read in bam_in.fetch():
        n_input += 1
        # Check if the read's reference name is in the set of contigs to include
        if read.reference_name in contigs_to_include:
            bam_out.write(read)
            n_written += 1

    # Close the BAM files
    bam_in.close()
    bam_out.close()

    return n_input, n_written

def main():
    args = parse_args()

    n_input, n_written = split_bam_by_contig(args.input, args.output, args.contig_file)

    if args.counts:
        with open(args.counts, "w") as f:
            f.write(f"{n_input}\t{n_written}\n")

if __name__ == "__main__":
    main()

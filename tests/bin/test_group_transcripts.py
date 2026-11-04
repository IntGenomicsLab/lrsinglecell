#!/usr/bin/env python3

""" Unit test for bin/group_transcripts.sh on three small FASTA/GTF fixtures.

    Not run by nf-test (which only collects *.nf.test); .github/workflows/bin-tests.yml
    runs it on pull requests. It needs only python3, bash and awk. The module's own
    container (ubuntu:20.04, whose awk is mawk) has no python, so the script is run in it
    through GROUP_TRANSCRIPTS_WRAPPER, a command prefix that must make $TMPDIR and the
    repository visible inside the container:

        python3 -m unittest -v tests/bin/test_group_transcripts.py

        GROUP_TRANSCRIPTS_WRAPPER="apptainer exec -B <repo> -B $TMPDIR <ubuntu:20.04 image>" \
            python3 -m unittest -v tests/bin/test_group_transcripts.py

    The fixtures cover the three header/annotation styles the script has to handle:
    GENCODE `|`-delimited headers with a GENCODE GTF, a versioned Ensembl / Cell Ranger
    GTF (`gene_version` before `transcript_id`, ids without a `.N` suffix), and Ensembl
    cDNA headers that contain spaces.
"""

import os
import shlex
import subprocess
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin", "group_transcripts.sh")
WRAPPER = shlex.split(os.environ.get("GROUP_TRANSCRIPTS_WRAPPER", ""))


def gtf_line(seqname, feature, attributes):
    return "\t".join([seqname, "HAVANA", feature, "100", "200", ".", "+", ".", attributes])


def gencode_tx(seqname, gene, tx):
    return [
        gtf_line(seqname, "transcript", f'gene_id "{gene}"; transcript_id "{tx}"; gene_type "lncRNA"; level 2;'),
        gtf_line(seqname, "exon", f'gene_id "{gene}"; transcript_id "{tx}"; gene_type "lncRNA"; exon_number 1; level 2;'),
    ]


def ensembl_tx(seqname, gene, tx, tx_version):
    attrs = (f'gene_id "{gene}"; gene_version "5"; transcript_id "{tx}"; transcript_version "{tx_version}"; '
             f'gene_type "protein_coding"; gene_name "G";')
    return [gtf_line(seqname, "transcript", attrs), gtf_line(seqname, "exon", attrs + " exon_number 1;")]


# GENCODE: `|`-delimited headers, versioned transcript_id as the second attribute, a chrY
# PAR copy whose id extends the chrX one, and a scaffold whose name contains a dot.
GENCODE_GTF = [
    "##description: evidence-based annotation of the human genome (GRCh38), version 45",
    "##provider: GENCODE",
    gtf_line("chr1", "gene", 'gene_id "ENSG00000290825.1"; gene_type "lncRNA"; gene_name "DDX11L2"; level 2;'),
    *gencode_tx("chr1", "ENSG00000290825.1", "ENST00000456328.2"),
    *gencode_tx("chr2", "ENSG00000227232.5", "ENST00000488147.1"),
    *gencode_tx("chrX", "ENSG00000182378.15", "ENST00000381192.10"),
    *gencode_tx("chrY", "ENSG00000182378.15_PAR_Y", "ENST00000381192.10_PAR_Y"),
    *gencode_tx("KI270734.1", "ENSG00000277196.4", "ENST00000615165.1"),
    *gencode_tx("chrM", "ENSG00000198888.2", "ENST00000361390.2"),
]
GENCODE_HEADERS = [
    ">ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|",
    ">ENST00000488147.1|ENSG00000227232.5|OTTHUMG00000000958.1|OTTHUMT00000002839.1|WASH7P-201|WASH7P|1351|unprocessed_pseudogene|",
    ">ENST00000381192.10|ENSG00000182378.15|OTTHUMG00000021073.6|-|PLCXD1-201|PLCXD1|2476|protein_coding|",
    ">ENST00000381192.10_PAR_Y|ENSG00000182378.15_PAR_Y|OTTHUMG00000021073.6|-|PLCXD1-201|PLCXD1|2476|protein_coding|",
    ">ENST00000615165.1|ENSG00000277196.4|-|-|AC007325.2-201|AC007325.2|1238|protein_coding|",
    ">ENST00000361390.2|ENSG00000198888.2|-|-|MT-ND1-201|MT-ND1|956|protein_coding|",
    ">ENST00000999999.1|ENSG00000999999.1|-|-|NOVEL-201|NOVEL|100|lncRNA|",
]
GENCODE_EXPECTED = {
    "chr1": [GENCODE_HEADERS[0][1:]],
    "chr2": [GENCODE_HEADERS[1][1:]],
    "chrX": [GENCODE_HEADERS[2][1:]],
    "chrY": [GENCODE_HEADERS[3][1:]],
    "KI270734.1": [GENCODE_HEADERS[4][1:]],
    "chrM": [GENCODE_HEADERS[5][1:]],
    "chr_Unknown": [GENCODE_HEADERS[6][1:]],
}

# Versioned Ensembl / Cell Ranger GTF: field 12 is the gene_version value, transcript_id is
# unversioned and the version sits in transcript_version. The FASTA ids carry the version.
ENSEMBL_GTF = [
    "#!genome-build GRCh38.p13",
    *ensembl_tx("1", "ENSG00000243485", "ENST00000473358", "1"),
    *ensembl_tx("1", "ENSG00000243485", "ENST00000469289", "1"),
    *ensembl_tx("10", "ENSG00000260370", "ENST00000562780", "1"),
    *ensembl_tx("MT", "ENSG00000198888", "ENST00000361390", "2"),
]
ENSEMBL_VERSIONED_HEADERS = [
    ">ENST00000473358.1",
    ">ENST00000562780.1",
    ">ENST00000469289.1",
    ">ENST00000361390.2",
]
ENSEMBL_VERSIONED_EXPECTED = {
    "1": ["ENST00000473358.1", "ENST00000469289.1"],
    "10": ["ENST00000562780.1"],
    "MT": ["ENST00000361390.2"],
}

# Ensembl cDNA: the header continues after a space, and the description contains `;`,
# which must not be taken as the id delimiter.
ENSEMBL_CDNA_HEADERS = [
    ">ENST00000473358.1 cdna chromosome:GRCh38:1:29554:31097:1 gene:ENSG00000243485.5 gene_biotype:lncRNA "
    "transcript_biotype:lncRNA gene_symbol:MIR1302-2HG description:MIR1302-2 host gene [Source:HGNC Symbol;Acc:HGNC:52482]",
    ">ENST00000562780.1 cdna chromosome:GRCh38:10:100:200:1 gene:ENSG00000260370.1 gene_biotype:lncRNA "
    "transcript_biotype:lncRNA",
    ">ENST00000361390.2 cdna chromosome:GRCh38:MT:3307:4262:1 gene:ENSG00000198888.2 gene_biotype:protein_coding "
    "transcript_biotype:protein_coding gene_symbol:MT-ND1 description:mitochondrially encoded NADH dehydrogenase 1 "
    "[Source:HGNC Symbol;Acc:HGNC:7455]",
]
ENSEMBL_CDNA_EXPECTED = {
    "1": ["ENST00000473358.1"],
    "10": ["ENST00000562780.1"],
    "MT": ["ENST00000361390.2"],
}


class GroupTranscriptsTest(unittest.TestCase):

    def run_script(self, headers, gtf, delimiter=None, expect_ok=True):
        """ Run the script in a fresh directory; return {chr: [lines]} of the files it wrote. """
        workdir = tempfile.mkdtemp(prefix="group_transcripts_", dir=os.environ.get("TMPDIR"))
        self.addCleanup(subprocess.run, ["rm", "-rf", workdir])
        with open(os.path.join(workdir, "in.fa"), "w") as f:
            for header in headers:
                f.write(header + "\nACGT\nACGT\n")
        with open(os.path.join(workdir, "in.gtf"), "w") as f:
            f.write("\n".join(gtf) + "\n")

        # cd inside the (possibly wrapped) shell, since a container does not inherit our cwd
        cmd = [*WRAPPER, "bash", "-c", 'cd "$0" && exec bash "$@"', workdir,
               os.path.abspath(SCRIPT), "-f", "in.fa", "-g", "in.gtf"]
        if delimiter is not None:
            cmd += ["-d", delimiter]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.stderr = result.stderr
        if expect_ok:
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertEqual(result.stdout, "", msg="the script should print nothing on stdout")
        else:
            self.assertNotEqual(result.returncode, 0)
            return None

        groups = {}
        for name in os.listdir(workdir):
            if name.endswith(".transcripts.txt"):
                with open(os.path.join(workdir, name)) as f:
                    groups[name[:-len(".transcripts.txt")]] = f.read().splitlines()
        return groups

    def test_gencode_auto_delimiter(self):
        self.assertEqual(self.run_script(GENCODE_HEADERS, GENCODE_GTF), GENCODE_EXPECTED)

    def test_gencode_explicit_delimiter(self):
        self.assertEqual(self.run_script(GENCODE_HEADERS, GENCODE_GTF, delimiter="|"), GENCODE_EXPECTED)

    def test_gencode_warns_about_unknown(self):
        self.run_script(GENCODE_HEADERS, GENCODE_GTF)
        self.assertIn("1 of 7 FASTA sequences matched no GTF transcript_id", self.stderr)

    def test_versioned_ensembl_gtf(self):
        self.assertEqual(self.run_script(ENSEMBL_VERSIONED_HEADERS, ENSEMBL_GTF), ENSEMBL_VERSIONED_EXPECTED)
        self.assertEqual(self.stderr, "")

    def test_versioned_gtf_unversioned_fasta(self):
        # GENCODE GTF ids carry a version, so an unversioned FASTA id must still find them
        headers = [">ENST00000456328", ">ENST00000361390"]
        self.assertEqual(self.run_script(headers, GENCODE_GTF),
                         {"chr1": ["ENST00000456328"], "chrM": ["ENST00000361390"]})

    def test_ensembl_cdna_headers_with_spaces(self):
        # The BAM reference name is the first token, and that is what has to be written
        self.assertEqual(self.run_script(ENSEMBL_CDNA_HEADERS, ENSEMBL_GTF), ENSEMBL_CDNA_EXPECTED)

    def test_ids_are_matched_exactly(self):
        # an unanchored grep for tx1 would also hit tx12 and tx1b
        gtf = [
            gtf_line("chrB", "transcript", 'gene_id "g12"; transcript_id "tx12";'),
            gtf_line("chrC", "transcript", 'gene_id "g1b"; transcript_id "tx1b";'),
            gtf_line("chrA", "transcript", 'gene_id "g1"; transcript_id "tx1";'),
        ]
        groups = self.run_script([">tx1", ">tx12", ">tx1b"], gtf)
        self.assertEqual(groups, {"chrA": ["tx1"], "chrB": ["tx12"], "chrC": ["tx1b"]})

    def test_no_match_at_all(self):
        groups = self.run_script([">unrelated_1", ">unrelated_2"], ENSEMBL_GTF)
        self.assertEqual(groups, {"chr_Unknown": ["unrelated_1", "unrelated_2"]})
        self.assertIn("none of the 2 FASTA sequences matched", self.stderr)

    def test_gtf_without_transcript_id_fails(self):
        gtf = [gtf_line("chr1", "gene", 'gene_id "g1"; gene_version "1";')]
        self.run_script(GENCODE_HEADERS, gtf, expect_ok=False)

    def test_missing_argument_fails(self):
        result = subprocess.run([*WRAPPER, "bash", os.path.abspath(SCRIPT), "-f", "/nonexistent.fa"],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()

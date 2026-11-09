# IntGenomicsLab/lrsinglecell: Tutorial

This tutorial takes you from a fresh machine to a Seurat object and per-cell DNA reads. It first runs the pipeline on the small public test dataset, so you can check your setup in a few minutes, and then shows how to run it on your own data.

The test dataset is a subsample of a 10X Multiome library from the RPMI-8402 cell line, sequenced on Oxford Nanopore. It has one cDNA (gene expression) FASTQ and one DNA (ATAC) FASTQ, aligned to chromosome 21 only. The test data is deliberately tiny (a few hundred reads), so the matrices it produces are not biologically meaningful. It exists to show you every step and every output file.

## 1. Install the prerequisites

You need:

1. **Nextflow ≥ 25.10.4**, which needs Java 17 or later. The simplest route is conda/mamba:

   ```bash
   mamba create -n nextflow -c bioconda -c conda-forge nextflow
   mamba activate nextflow
   nextflow -version
   ```

   or follow the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html).

2. **One container engine or package manager** that provides the tools. Pick whichever your system has; you pass it to the pipeline with `-profile`:

   | You have                                      | Use                                            |
   | --------------------------------------------- | ---------------------------------------------- |
   | Docker (laptop / workstation)                 | `-profile docker`                              |
   | Singularity or Apptainer (most HPC clusters)  | `-profile singularity` or `-profile apptainer` |
   | Neither (slower to set up, less reproducible) | `-profile conda`                               |

   You do **not** need to install minimap2, flexiplex, IsoQuant or any other tool yourself; every process runs in its own container.

On an HPC cluster, first check whether your institute has a profile in [nf-core/configs](https://github.com/nf-core/configs). If it does, add it, e.g. `-profile singularity,<institute>`, and the pipeline will submit its jobs to your scheduler.

## 2. Run the test profile

```bash
nextflow run IntGenomicsLab/lrsinglecell \
    -r v1.0.1 \
    -profile test_mixed,singularity \
    --outdir results_test
```

- `-r v1.0.1` pins a released version, so you get the same pipeline every time. Always pin a version for real analyses.
- `test_mixed` fills in all required parameters (samplesheet, chr21 reference, GTF, `--barcode_format 10X_multiome`, `--quantifier isoquant`) from the public [test-datasets](https://github.com/IntGenomicsLab/test-datasets/tree/main/lrsinglecell) repository. `test` (cDNA only) and `test_dna` (DNA only) are also available.
- The first run downloads the containers, which takes most of the time. On our cluster the whole run takes about 7 minutes.

A successful run ends with:

```text
-[IntGenomicsLab/lrsinglecell] Pipeline completed successfully-
```

If it fails, the error message names the failing process and its work directory. `.command.err` and `.command.log` in that directory contain the tool's own output.

## 3. Write your samplesheet

The test profile used [this samplesheet](https://raw.githubusercontent.com/IntGenomicsLab/test-datasets/main/lrsinglecell/samplesheet/samplesheet_test_mixed.csv):

```csv title="samplesheet.csv"
sample,fastq,cell_count,type
RPMI8402,https://raw.githubusercontent.com/IntGenomicsLab/test-datasets/main/lrsinglecell/fastq/sub_RPMI8402_cDNA.fastq.gz,1000,cdna
RPMI8402,https://raw.githubusercontent.com/IntGenomicsLab/test-datasets/main/lrsinglecell/fastq/sub_RPMI8402_DNA.fastq.gz,1000,dna
```

For your own data, write one row per FASTQ file:

| Column       | Required | What to put there                                                                                                                           |
| ------------ | -------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`     | yes      | Sample name without spaces. Give the cDNA and DNA libraries of one Multiome sample the same name; they are processed separately, by `type`. |
| `fastq`      | yes      | Absolute path (or URL) to a gzipped ONT FASTQ (`.fastq.gz` / `.fq.gz`). cDNA and DNA libraries are separate files.                          |
| `cell_count` | yes      | Expected number of cells. Only used by the BLAZE barcode caller, but must be filled in on every row.                                        |
| `type`       | no       | `cdna` (default) or `dna`.                                                                                                                  |

If one library was sequenced on several flow cells, add one row per FASTQ with the same `sample` and `type`; the reads are concatenated before processing. If your basecaller wrote many small FASTQ files per flow cell, concatenate them first (`cat *.fastq.gz > run1.fastq.gz`).

## 4. Run on your own data

You need a genome FASTA and (for cDNA) a matching GTF, e.g. from [GENCODE](https://www.gencodegenes.org/human/) or a 10X Cell Ranger reference.

10X Multiome (cDNA + DNA):

```bash
nextflow run IntGenomicsLab/lrsinglecell -r v1.0.1 \
    -profile singularity \
    --input samplesheet.csv \
    --genome_fasta /refs/GRCh38.primary_assembly.genome.fa \
    --gtf /refs/gencode.v45.annotation.gtf \
    --barcode_format 10X_multiome \
    --quantifier isoquant \
    --outdir results
```

The settings to change for your data:

| You have                        | Set                                                                                                                     |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| 10X 3' v3 / v4 or 5' v2 / v3    | `--barcode_format 10X_3v3`, `10X_3v4`, `10X_5v2` or `10X_5v3`                                                           |
| 10X Multiome                    | `--barcode_format 10X_multiome` (required for DNA rows unless you give a custom pattern)                                |
| Another DNA barcoding chemistry | `--custom_flexiplex_barcode_dna "<flexiplex pattern>" --dna_whitelist barcodes.txt` ([details](usage.md#dna-libraries)) |
| Transcript quantification only  | `--quantifier oarfish --transcript_fasta transcripts.fa`, or `--quantifier lrkallisto`                                  |
| Several quantifiers             | `--quantifier isoquant,oarfish,lrkallisto`                                                                              |
| DNA only                        | omit `--gtf` and `--quantifier`                                                                                         |

To keep the command short and reproducible, put the parameters in a YAML file and pass it with `-params-file params.yml`. Add `-resume` to restart a failed or interrupted run without redoing finished steps. For large PromethION runs, see the resource notes in the [README troubleshooting section](../README.md#troubleshooting).

All parameters are listed with `nextflow run IntGenomicsLab/lrsinglecell -r v1.0.1 --help_full` and documented in the [usage docs](usage.md).

## 5. Find your way around the results

This is the output of the test run above (`results_test/`), trimmed to the files you will use most:

```text
results_test/
├── RPMI8402/
│   ├── cdna/
│   │   ├── flexiplex/                      barcode counts and the called-cell barcode list
│   │   ├── genome/
│   │   │   ├── bam/dedup/                  RPMI8402.cdna.genome.dedup.bam(.bai): genome alignment, UMI-deduplicated
│   │   │   ├── isoquant/RPMI8402.cdna/RPMI8402.cdna/
│   │   │   │   ├── RPMI8402.cdna.gene_grouped_tag_CB_counts.{matrix.mtx,barcodes.tsv,features.tsv}
│   │   │   │   └── RPMI8402.cdna.transcript_grouped_tag_CB_counts.{matrix.mtx,barcodes.tsv,features.tsv}
│   │   │   ├── isoquant_all_droplets/      the same, for every droplet (for ambient-RNA tools)
│   │   │   └── qc/seurat_isoquant/gene/    RPMI8402_seurat.rds, QC plots and statistics
│   │   └── qc/nanoplot/                    read-level QC before/after trimming
│   └── dna/
│       ├── flexiplex/                      barcode counts and the called-cell barcode list
│       ├── bam/dedup/RPMI8402.bam(.bai)    genome alignment, duplicates flagged
│       └── qc/dedup/                       duplication metrics
├── batch_qcs/multiqc/final_qc/multiqc_report.html   start here: all QC in one report
├── references/                             minimap2 indexes, reusable with -resume
└── pipeline_info/                          run reports, execution trace, software versions
```

Start with `multiqc_report.html`. [output.md](output.md) describes every file.

Every alignment carries its cell barcode in BAM tags:

| Tag  | Meaning                                                                              |
| ---- | ------------------------------------------------------------------------------------ |
| `CB` | Cell barcode corrected against the whitelist; `-` if the read matched no called cell |
| `CR` | Raw barcode as read from the sequence                                                |
| `XB` | `CB` when there is one, otherwise `CR`: the droplet the read came from               |
| `UB` | UMI (cDNA only)                                                                      |
| `GX` | Gene the read was assigned to (cDNA only)                                            |

## 6. Load the results

### Gene expression in Seurat (R)

The pipeline already saves a Seurat object per sample and quantifier:

```r
library(Seurat)

obj = readRDS("results_test/RPMI8402/cdna/genome/qc/seurat_isoquant/gene/RPMI8402_seurat.rds")
obj
```

Or build one yourself from the matrix, e.g. to combine samples:

```r
library(Seurat)

mtx_dir = "results_test/RPMI8402/cdna/genome/isoquant/RPMI8402.cdna/RPMI8402.cdna"
prefix = file.path(mtx_dir, "RPMI8402.cdna.gene_grouped_tag_CB_counts")

counts = ReadMtx(
  mtx = paste0(prefix, ".matrix.mtx"),
  cells = paste0(prefix, ".barcodes.tsv"),
  features = paste0(prefix, ".features.tsv"),
  feature.column = 1
)
obj = CreateSeuratObject(counts, project = "RPMI8402")
```

From here, follow any standard Seurat workflow. The transcript-level matrix (`transcript_grouped_tag_CB_counts`) loads the same way and can be added as a second assay for isoform analysis. The lr-kallisto `gene/` directory is in Cell Ranger format and loads directly with `Read10X()`.

In Python, `scanpy.read_mtx()` reads the same files (transpose to cells × genes).

### Per-cell DNA reads

The DNA BAM holds all cells of a sample, each read tagged with its cell barcode. Duplicates are flagged (`0x400`), not removed. To get the deduplicated reads of one cell:

```bash
samtools view -b -F 0x400 -d CB:GAGTAACTCTCAATAC \
    results_test/RPMI8402/dna/bam/dedup/RPMI8402.bam > cell_GAGTAACTCTCAATAC.bam
```

or to split the whole sample into one BAM per called cell. The called-cell barcodes (with their read counts) are listed in `flexiplex/<sample>.dna_known_barcodes.txt`:

```bash
mkdir -p cells
while read -r barcode count; do
    samtools view -b -F 0x400 -d CB:"$barcode" \
        results_test/RPMI8402/dna/bam/dedup/RPMI8402.bam > cells/"$barcode".bam
done < results_test/RPMI8402/dna/flexiplex/RPMI8402.dna_known_barcodes.txt
```

For a 10X Multiome sample, the DNA and cDNA barcodes of the same cell are different sequences, taken from the ATAC and gene-expression whitelists. The two whitelists bundled in `assets/whitelist/` (`cellranger_arc_atac.737K-arc-v1.txt.gz` and `cellranger_arc_rna.737K-arc-v1.txt.gz`) are paired line by line: the barcode on line _n_ of one belongs to the same gel bead as line _n_ of the other.

## Getting help

- Parameter reference: [usage.md](usage.md) and `--help`
- Output reference: [output.md](output.md)
- Questions, bugs and feature requests: [GitHub issues](https://github.com/IntGenomicsLab/lrsinglecell/issues)

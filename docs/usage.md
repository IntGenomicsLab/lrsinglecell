# IntGenomicsLab/lrsinglecell: Usage

> _The full list of parameters, with their defaults, is in the pipeline schema ([`nextflow_schema.json`](../nextflow_schema.json)) and is printed by `nextflow run IntGenomicsLab/lrsinglecell --help` (add `--help_full` for every parameter, `--show_hidden` for hidden ones). This page explains how to use the parameters that matter most._

A step-by-step walkthrough, from installation to loading the results in R, is in the [tutorial](tutorial.md). The files the pipeline writes are described in [output.md](output.md).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with a header row, three required columns and an optional fourth, as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

The example `samplesheet.csv` below describes four 10X Multiome samples, each with one cDNA and one DNA FASTQ file, and sample-specific cell counts. The cDNA and DNA libraries of one sample are separate FASTQ files.

```csv title="samplesheet.csv"
sample,fastq,cell_count,type
CONTROL_REP1,/data/AEG588A1_cDNA.fastq.gz,5000,cdna
CONTROL_REP2,/data/AEG588A2_cDNA.fastq.gz,6000,cdna
TREATMENT_REP1,/data/AEG588A4_cDNA.fastq.gz,5500,cdna
TREATMENT_REP2,/data/AEG588A5_cDNA.fastq.gz,6000,cdna
CONTROL_REP1,/data/AEG588A1_DNA.fastq.gz,5000,dna
CONTROL_REP2,/data/AEG588A2_DNA.fastq.gz,6000,dna
TREATMENT_REP1,/data/AEG588A4_DNA.fastq.gz,5500,dna
TREATMENT_REP2,/data/AEG588A5_DNA.fastq.gz,6000,dna
```

| Column       | Description                                                                                                                                                                                  |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`     | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. It must not contain spaces.                                                    |
| `fastq`      | Full path (or URL) to a FastQ file for Oxford Nanopore. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". Each file can appear only once in the samplesheet.            |
| `cell_count` | Expected number of cells/nuclei, an integer greater than 0. Only BLAZE (`--demux_tool_cdna blaze`) uses it, as `--expect-cells`. It is required on DNA rows too, but DNA rows do not use it. |
| `type`       | An optional column specifying whether the library is DNA (`dna`) or cDNA (`cdna`). If omitted, the default `cdna` is used.                                                                   |

Rows with the same `sample` and `type` must use the same `cell_count`; the run stops with an error otherwise.

Note: DNA samples are only compatible with `flexiplex` demultiplexing. They also need `--genome_fasta` and a known barcode layout; see [DNA libraries](#dna-libraries).

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline, and the samplesheets used by the test profiles (with real, downloadable FASTQ URLs) are in [IntGenomicsLab/test-datasets](https://github.com/IntGenomicsLab/test-datasets/tree/main/lrsinglecell/samplesheet).

### DNA libraries

DNA rows (`type` = `dna`) are demultiplexed with flexiplex and need a barcode layout plus a matching whitelist. There are two ways to give them:

1. **10X Multiome (ATAC side)**: `--barcode_format 10X_multiome`. The pipeline then uses the bundled ATAC whitelist `assets/whitelist/cellranger_arc_atac.737K-arc-v1.txt.gz` for DNA, and the RNA whitelist `assets/whitelist/cellranger_arc_rna.737K-arc-v1.txt.gz` for any cDNA rows in the same run. The flexiplex pattern used is

   ```text
   -x ACCGAGATCTACAC -b ???????????????? -x CGCGTCTGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -f 8
   ```

   (end of the P5 adapter, 16 bp cell barcode, then the Nextera read 1 + Tn5 mosaic-end sequence; DNA libraries carry no UMI).

2. **Any other barcoding chemistry**: give your own flexiplex pattern with `--custom_flexiplex_barcode_dna` and a plain-text or gzipped whitelist (one barcode per line) with `--dna_whitelist`. For example:

   ```bash
   --custom_flexiplex_barcode_dna "-x ACCGAGATCTACAC -b ???????????????? -x CGCGTCTGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -f 8" \
   --dna_whitelist /path/to/barcodes.txt.gz
   ```

   See the [flexiplex documentation](https://davidsongroup.github.io/flexiplex/) for the pattern syntax.

The pipeline checks this at start-up: a samplesheet with DNA rows but no known layout, or a custom pattern without `--dna_whitelist`, stops with an error before any job is submitted. DNA-only runs do not need `--gtf` or `--quantifier`.

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads of rows with the same `sample` and `type` before performing any downstream analysis. Below is an example for the same sample sequenced across replicates 1 and 4 (`REP1` and `REP4` respectively):

```csv title="samplesheet.csv"
sample,fastq,cell_count,type
CONTROL_REP1,AEG588A1_S1.fastq.gz,5000,cdna
CONTROL_REP1,AEG588A1_S2.fastq.gz,5000,cdna
CONTROL_REP2,AEG588A2_S1.fastq.gz,2000,cdna
CONTROL_REP3,AEG588A3_S1.fastq.gz,7500,cdna
CONTROL_REP4,AEG588A4_S1.fastq.gz,9000,cdna
CONTROL_REP4,AEG588A4_S2.fastq.gz,9000,cdna
CONTROL_REP4,AEG588A4_S3.fastq.gz,9000,cdna
```

### Dorado duplex reads

`dorado duplex` writes three kinds of reads to the same output, told apart by the `dx` tag:

| `dx`      | read                                                                    |
| --------- | ----------------------------------------------------------------------- |
| `dx:i:1`  | duplex consensus read, built from the template and complement strand    |
| `dx:i:0`  | simplex read without a duplex partner                                   |
| `dx:i:-1` | simplex template or complement read that also went into a `dx:i:1` read |

Without filtering, each duplex molecule enters the pipeline three times: the consensus read and both of its parents. Deduplication collapses most of these, but with `--drop_duplex_parents` the `dx:i:-1` reads are removed before QC and demultiplexing, so each molecule is used once, as its more accurate consensus read, with its barcode read from both strands. The option is off by default.

```bash
nextflow run IntGenomicsLab/lrsinglecell ... --drop_duplex_parents
```

The FASTQ headers must carry dorado's tags. **`dorado duplex --emit-fastq` does not write the `dx` tag** (checked with dorado 1.2.0), so convert dorado's duplex BAM instead, keeping the tags with `-T '*'`:

```bash
samtools fastq -T '*' calls.bam | gzip > reads.fastq.gz
```

A plain `samtools fastq` without `-T '*'` drops the tags too.

With `--drop_duplex_parents` set, the pipeline looks at the first 1000 reads of each sample:

- **`dx` tags found**: the `dx:i:-1` reads are dropped and the counts are written to `<sample>/<type>/qc/duplex_filter/` (see the [output docs](output.md#duplex-filter)). If the file has `dx:i:-1` reads but no `dx:i:1` reads, the consensus reads were most likely removed upstream; all reads are then kept, with a warning, as dropping the parents would lose these molecules.
- **No `dx` tags, but duplex consensus reads (read names `<id1>;<id2>`)**: the run stops with an error. This is what `--emit-fastq` output looks like, and the parents cannot be identified in it. Convert from BAM as above, or run without `--drop_duplex_parents`.
- **No `dx` tags and no duplex reads** (simplex data): the FASTQ is passed through unchanged, without reading it in full, and the pipeline prints a warning.

Because only the first reads are checked, do not mix simplex and duplex FASTQs within one sample. A duplex FASTQ is rewritten in a single pass with 8 CPUs, which recompress the kept reads with `pigz`. This took about 8 s per GB of uncompressed FASTQ in our tests, so very large inputs (hundreds of GB of `fastq.gz`) take one to two hours.

## Running the pipeline

The typical command for running the pipeline on a 10X Multiome samplesheet (cDNA + DNA) is as follows:

```bash
nextflow run IntGenomicsLab/lrsinglecell -r <version> \
  --input ./samplesheet.csv \
  --outdir ./results \
  --genome_fasta /path/to/genome.fa \
  --gtf /path/to/file.gtf \
  --quantifier isoquant,oarfish \
  --demux_tool_cdna flexiplex \
  --demux_tool_dna flexiplex \
  --barcode_format 10X_multiome \
  -profile <docker/singularity/apptainer/institute>
```

Replace `<version>` with a release tag, e.g. `-r v1.1.0` (see [Reproducibility](#reproducibility)).

For a cDNA-only run, set `--barcode_format` to the 10X chemistry of the library (`10X_3v3`, `10X_3v4`, `10X_5v2` or `10X_5v3`); this also selects the matching bundled cDNA whitelist, which you can override with `--cdna_whitelist`. See [Barcode formats](#barcode-formats).

Which references are needed depends on the data types and quantifiers:

- `--genome_fasta` is required for DNA samples, for the `isoquant` and `lrkallisto` quantifiers, and for `oarfish` when `--transcript_fasta` is not given.
- `--gtf` is required whenever the samplesheet contains `cdna` rows, and can be omitted for a DNA-only run.
- `--transcript_fasta` is optional and only used by `oarfish`. Without it, the transcript sequences are extracted from `--genome_fasta` and `--gtf` with `gffread -w` (the same transcriptome lr-kallisto indexes) and published to `references/transcriptome/transcripts.fasta`. If you give one, e.g. GENCODE `transcripts.fa`, its sequence ids must match the `transcript_id`s of the GTF (after [`--fasta_delimiter`](#transcriptome-deduplication-and---fasta_delimiter)). lr-kallisto always uses the `gffread` transcriptome.

For the `--quantifier` parameter, give one quantifier (`isoquant`, `oarfish` or `lrkallisto`) or several as a comma-separated list without spaces, e.g. `isoquant,oarfish`.

### Quantifiers

Three quantifiers are available through `--quantifier`, and any combination of them can be requested as a comma-delimited string:

| Quantifier   | Reference required                                            | Levels            | Notes                                                                  |
| ------------ | ------------------------------------------------------------- | ----------------- | ---------------------------------------------------------------------- |
| `isoquant`   | `genome_fasta` + `gtf`                                        | gene + transcript | Quantifies from a genome alignment                                     |
| `oarfish`    | `gtf`, plus `transcript_fasta` or `genome_fasta` to derive it | transcript        | Quantifies from a transcriptome alignment; always deduplicates         |
| `lrkallisto` | `genome_fasta` + `gtf`                                        | gene + transcript | Pseudoaligns straight from the FASTQ; quantification uses no alignment |

`--quantifier` only applies to `cdna` samples. It is required whenever the samplesheet contains at least one `cdna` row, and can be omitted entirely for a DNA-only run: DNA samples are aligned and deduplicated but never quantified. If a DNA-only samplesheet is given without `--quantifier`, the pipeline runs the DNA path only; if a `cdna` row is present without `--quantifier`, the run stops with an error.

The cDNA genome alignment is made when `isoquant` or `lrkallisto` is requested, and the transcriptome alignment when `oarfish` is. An `oarfish`-only run therefore has no genome BAM.

#### lr-kallisto

[lr-kallisto](https://kallisto.readthedocs.io/en/latest/lr/pseudoalignment.html) is the long-read mode of kallisto (`kallisto --long`). Unlike the other two quantifiers, it does not align reads: it pseudoaligns them against a k=63 index and quantifies with a long-read expectation-maximisation step. That makes the quantification itself considerably faster.

The cDNA genome alignment is still produced when `lrkallisto` is requested, so a BAM and its alignment-derived QC are available under `<sample>/cdna/genome/` even on an lr-kallisto-only run. lr-kallisto does not read it. If you want lr-kallisto's speed and do not need the deduplicated BAM, `--skip_cdna_dedup` drops the expensive `umi_tools` chain and publishes the mapped-only BAM instead.

Because lr-kallisto locates the cell barcode and UMI positionally rather than by name, the pipeline moves them out of the flexiplex read name into a synthetic barcode read before handing the FASTQs to `kallisto bus`. This reuses the barcodes that flexiplex already corrected against the whitelist, so cell barcodes are directly comparable with the other quantifiers. UMI deduplication is performed by `bustools` while counting, so `--dedup_tool`, `--skip_cdna_dedup` and `--skip_dedup` have no effect on the lr-kallisto matrices. The published `gene/` matrix is an integer UMI count matrix (`bustools count --genecounts`) and can be given directly to CellBender or any other tool that expects counts; the `transcript/` matrix and the gene abundances under `quant/` are `kallisto quant-tcc` EM estimates and are fractional. See [`docs/output.md`](output.md#lr-kallisto) for the difference.

This means lr-kallisto has some extra requirements:

- `--genome_fasta` and `--gtf` must both be provided. The transcript sequences are extracted from them with `gffread`, the transcript-to-gene mapping is derived from the same GTF, and the genome is used as the kallisto d-list so that reads originating outside the transcriptome are discarded rather than misassigned. Deriving both from one GTF is what guarantees the index sequence names and the mapping agree.
- `--demux_tool_cdna` must be `flexiplex`; the `blaze` path keeps barcodes in a separate per-read table rather than in the read name.
- `--custom_flexiplex_barcode_cdna` cannot be used, because the barcode and UMI widths are needed in order to locate them. Use one of the `--barcode_format` presets instead.

lr-kallisto is run twice: once on the called cells (`lrkallisto/`) and once on every droplet (`lrkallisto_all_droplets/`), for ambient-RNA and empty-droplet tools. Two parameters control the second pass:

- `--lrkallisto_all_min_reads` (default `100`) keeps only the barcodes flexiplex saw at least that many times in the all-droplet pass. Without it, sequencing errors in the uncorrected below-knee barcodes create millions of near-singleton columns, and the EM time grows with them. Called cells sit far above the threshold. Set it to `0` to keep every barcode.
- `--skip_lrkallisto_all_droplets` skips the all-droplet pass entirely, halving the lr-kallisto counting work. Nothing else in the pipeline uses its output.

The k-mer length, platform and d-list can be tuned with `--kallisto_kmer_size` (default 63), `--kallisto_platform` (`ONT` or `PacBio`, default `ONT`) and `--kallisto_dlist` (default `true`). `--kallisto_threshold` (default 0.8) is also accepted, but it is inert; see the last warning in this section.

> [!WARNING]
> To turn the d-list off, write `--kallisto_dlist=false` with an equals sign. `params.kallisto_dlist` defaults to a boolean, so Nextflow treats the space-separated `--kallisto_dlist false` as a bare flag and drops the `false` as a positional argument: the d-list stays **on** and the only sign is a `positional argument` warning in the log. The same applies to any other boolean parameter you want to set to `false` on the command line.

> [!IMPORTANT]
> Check `p_pseudoaligned` in `bus/run_info.json` after your first run. `--kallisto_dlist` costs real sensitivity; `--kallisto_kmer_size` is already at its best value:
>
> |                                     | `--kallisto_dlist=false` | `--kallisto_dlist=true` (default) |
> | ----------------------------------- | ------------------------ | --------------------------------- |
> | `--kallisto_kmer_size 31`           | 31.5%                    | 18.8%                             |
> | `--kallisto_kmer_size 63` (default) | **42.8%**                | 30.1%                             |
>
> Measured on 9,805,796 flexiplex-demultiplexed ONT reads against GRCh38 (CellRanger `refdata-gex-GRCh38-2020-A`). On the same reads minimap2 places 95% and isoquant then assigns 45.8% of them to a gene, so 42.8% is in the same range as the alignment-based path rather than an order of magnitude below it.
>
> - **k-mer length.** Pseudoalignment needs _exact_ k-mer matches, so k interacts with per-base error rate, but **k=63 beats k=31 on real data at both d-list settings**, and it wins by most on the reads that matter: of the reads isoquant assigns confidently to a single gene, k=63 recovers 88.6% against k=31's 58.3% (d-list off). Do not lower k to chase a mapping rate. (An earlier version of this table reported the opposite ordering from the chr21 `test_lrkallisto` profile. That reference is too small for the comparison to carry: every rate there was under 21%, and the k ordering inverts on a whole transcriptome.)
> - **d-list.** Using the genome as a d-list discards reads containing k-mers that are in the genome but not the transcriptome, which is what stops intronic and genomic reads being force-assigned to a transcript. This is far more punitive for long reads than for short ones: a long read only has to contain a single genome-only k-mer to be rejected. Turning it off recovered 1,253,671 extra reads. Classifying each of them by what isoquant made of the same read: 59% are reads isoquant also assigns to a gene, 22% are intronic reads isoquant places inside a gene but does not count, 14% are reads minimap2 never placed at all, and only **3% are intergenic**, the one class that is unambiguously spurious. On this data the d-list costs far more sensitivity than the specificity it buys, so `--kallisto_dlist=false` is worth trying for ONT cDNA. The pipeline default is still `true`, which is what lr-kallisto upstream recommends; this was one flow cell of one chemistry, so it is a reason to test the setting on your own data rather than to assume the default is wrong.
> - **read length is not the limiting factor.** Only 1.1% of these reads were shorter than k=63, so reads too short to carry a single k-mer explain almost none of the loss.
>
> `--kallisto_threshold` will not recover anything: it is **inert**.
>
> The defaults are deliberately the upstream-recommended, higher-specificity ones. If your mapping rate is far below the rate you get from minimap2 on the same data, compare against another quantifier before trusting the counts.

> [!WARNING]
> `--kallisto_threshold` currently has no effect. `kallisto bus --long` accepts `--threshold` and silently ignores it: sweeping it over 0.0, 0.1, 0.2, 0.4, 0.6, 0.8 and 0.9 on 1,000,000 reads returned the identical set of 427,251 pseudoaligned reads **and** the identical resolved transcript set for every one of them; only the equivalence-class numbering moved. `--threshold` appears in no other kallisto subcommand, and the pipeline quantifies with `kallisto quant-tcc --long`, so the parameter cannot take effect anywhere on this path. It is retained only so that existing configs keep parsing; do not use it to tune sensitivity.

Reads are pseudoaligned as non-strand-specific (`kallisto bus --unstranded`), and `--stranded` does not apply to this path. Whether the reads flexiplex writes out are in a consistent orientation, which a stranded mode would need, is being investigated in [#46](https://github.com/IntGenomicsLab/lrsinglecell/issues/46). If you do want to force a strand, override the `KALLISTO_BUS` arguments with a custom config. The process runs as `...:QUANTIFY_SCRNA_LRKALLISTO:LRKALLISTO_COUNT:KALLISTO_BUS` (called cells) and `...:QUANTIFY_SCRNA_LRKALLISTO:LRKALLISTO_COUNT_ALL:KALLISTO_BUS` (all droplets); this selector matches both:

```groovy title="custom.config"
process {
    withName: '.*:LRKALLISTO_COUNT(_ALL)?:KALLISTO_BUS' {
        ext.args = { "--threshold ${params.kallisto_threshold} --fr-stranded" }
    }
}
```

> [!NOTE]
> Building the index is cheaper than you might expect for a k=63 index with a whole-genome d-list. Measured on GRCh38 (CellRanger `refdata-gex-GRCh38-2020-A`, 199,138 transcripts, 345 MB of transcript sequence) with the 3.15 Gb genome as the d-list: **463 s wall on 24 cpus and 3.6 GB peak RSS**, producing a 222 MB index. The `process_high` default is comfortably sufficient; there is no need to route `KALLISTO_INDEX` to a high-memory queue. If you use a much larger reference, check `peak_rss` for that process in the execution trace before assuming these numbers carry over.

### Barcode formats

The pipeline supports barcode identification and extraction through both `flexiplex` (default) and `blaze`, set with `--demux_tool_cdna`. DNA is always demultiplexed with `flexiplex` (`--demux_tool_dna flexiplex`). `--barcode_format` is required and selects, per data type, the bundled whitelist and the flexiplex pattern:

| `--barcode_format` | Bundled cDNA whitelist (`assets/whitelist/`) | Bundled DNA whitelist                    | flexiplex cDNA pattern                                                                   | Barcode / UMI (bp)          | Demultiplexers                                               |
| ------------------ | -------------------------------------------- | ---------------------------------------- | ---------------------------------------------------------------------------------------- | --------------------------- | ------------------------------------------------------------ |
| `10X_3v3`          | `3M-february-2018.txt.gz`                    | none                                     | `-d 10x3v3`                                                                              | 16 / 12                     | cDNA: flexiplex, BLAZE (`--kit-version 3v3`). No DNA         |
| `10X_3v4`          | `3M-3pgex-may-2023_TRU.txt.gz`               | none                                     | `-d 10x3v3` (same geometry as 3' v3)                                                     | 16 / 12                     | cDNA: flexiplex, BLAZE (`--kit-version 3v4`). No DNA         |
| `10X_5v2`          | `737K-august-2016.txt.gz`                    | none                                     | `-d 10x5v2`                                                                              | 16 / 10                     | cDNA: flexiplex, BLAZE (`--kit-version 5v2`). No DNA         |
| `10X_5v3`          | `3M-5pgex-jan-2023.txt.gz`                   | none                                     | `-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTCTTATATGGG -f 8 -e 2` | 16 / 10                     | cDNA: flexiplex, BLAZE (`--kit-version 5v3`). No DNA         |
| `10X_multiome`     | `cellranger_arc_rna.737K-arc-v1.txt.gz`      | `cellranger_arc_atac.737K-arc-v1.txt.gz` | `-d 10x3v3`; DNA: see [DNA libraries](#dna-libraries)                                    | cDNA 16 / 12; DNA 16 / none | cDNA: flexiplex, BLAZE (`--kit-version 3v3`). DNA: flexiplex |

Notes:

- `--cdna_whitelist` and `--dna_whitelist` replace the bundled whitelist of their data type. A whitelist is a plain-text or gzipped file with one barcode per line.
- `--custom_flexiplex_barcode_cdna` and `--custom_flexiplex_barcode_dna` replace the flexiplex pattern of their data type; `--barcode_format` still selects the bundled whitelist. A custom pattern is a string of flexiplex options, e.g. `"-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTCTTATATGGG -f 8 -e 2"`; see the [flexiplex documentation](https://davidsongroup.github.io/flexiplex/). Custom patterns only apply to flexiplex, and a custom cDNA pattern cannot be combined with `lrkallisto`.
- DNA rows need `--barcode_format 10X_multiome`, or a custom DNA pattern together with `--dna_whitelist`.
- With BLAZE, `10X_multiome` cDNA is processed as 3' v3. The barcode and UMI lengths in the table are also the ones the BLAZE path and lr-kallisto use to cut out the barcode and UMI.

### Read trimming and filtering

Reads are filtered with [chopper](https://github.com/wdecoster/chopper) before barcode calling:

- `--min_q_score` (default `10`): minimum average read quality (`chopper --quality`).
- `--min_length` (default `1`): minimum read length (`chopper --minlength`).
- `--skip_trimming`: skip chopper altogether; the post-trim read QC is then not produced.

With `--drop_duplex_parents` (off by default), the simplex parents of dorado duplex reads are removed before the pre-trim read QC; see [Dorado duplex reads](#dorado-duplex-reads).

### Strandedness

Everything runs unstranded by default: `--stranded` is unset, so minimap2 is run with `-un` and IsoQuant with `--stranded none`. `--stranded forward` or `--stranded reverse` switch both to the given strand (`-uf`/`-ub`, `--stranded forward/reverse`); lr-kallisto is always unstranded (see above). Which setting is right for each 10X kit has not been established yet; see [#46](https://github.com/IntGenomicsLab/lrsinglecell/issues/46).

### Deduplication

cDNA and DNA are deduplicated differently, and each arm has its own switch:

- **cDNA** reads are UMI-deduplicated with UMI-tools (`--dedup_tool umitools`, the only option; the Picard alternative was removed in v1.2.0). The transcriptome alignment used by oarfish is always deduplicated. Deduplication of the genome alignment (used by IsoQuant) can be turned off with `--skip_cdna_dedup`, in which case the mapped-only BAM is published and quantified instead.
- **DNA** reads carry no UMI and are duplicate-marked by `bin/mark_dna_duplicates.py`. `--skip_dedup` turns this off; it has no effect on cDNA. On the contigs in `--dna_exact_only_contigs` (comma-separated, default `chrM,MT`, matched exactly) only the exact duplicate rules are applied, because the mitochondrial genome's dense coverage makes the position-tolerant rules merge distinct molecules. If your reference names the mitochondrial genome differently, give its name here; the pipeline logs a warning when none of the listed contigs is in the reference (`.fai`). An empty value (`dna_exact_only_contigs: ""` in a [parameter file](#parameter-files)) applies the tolerant rules everywhere. See [DNA duplicate marking](output.md#dna-duplicate-marking) for the rules.

By default UMI-tools groups UMIs by feature rather than by alignment position: `--dedup_per_gene` (genome alignments) and `--dedup_per_contig` (transcriptome alignments), both `true`. This matters because `umi_tools` corrects its bundling position for soft clipping, and long reads carry variable amounts of untrimmed adapter and polyA, so reads from one PCR family end up in different position groups even when they align at an identical coordinate. Position-based grouping consequently leaves a large fraction of duplicates in place, and absolute molecule counts come out inflated.

- `--dedup_per_gene` tags every alignment with the gene whose body it overlaps most (`GX`/`GN`/`GS` tags, written by `TAG_GENES`) and groups on that gene. It requires `--gtf`. Reads whose gene call is ambiguous, or that overlap no gene at all, are deduplicated by position instead and merged back in, so no reads are dropped. The `GS` distribution per contig, which includes the reads without a gene, is published under `qc/gene_assignment/`, next to `*.gene_tag_summary.tsv`, which holds the per-sample `GS=none` fraction (also shown in MultiQC as "% reads without gene"; see [Gene assignment](output.md#gene-assignment)). The genome fasta and the gtf must use the same contig names: `TAG_GENES` fails if the BAM header shares no contig with the gtf (e.g. `chr1` vs `1`), and warns if fewer than half of the gtf's gene-bearing contigs are in the BAM.
- `--dedup_per_contig` applies to transcriptome alignments only, where a contig is a transcript. It needs no gene tag and no fallback.
- Set either to `false` to restore position-based grouping. Note that these change absolute UMI counts substantially. Relative expression is essentially unaffected, but anything quoting UMIs per cell, library complexity or saturation will differ from earlier releases.

### Analysis options

- `--retain_introns` (default `true`): count intronic and otherwise inconsistent reads in the IsoQuant matrices. With `true`, IsoQuant runs with `--gene_quantification all --transcript_quantification all`; with `--retain_introns=false` it uses `--gene_quantification unique_inconsistent --transcript_quantification unique_only`. Single-nucleus data usually needs the default.
- `--save_genome_secondary_alignment` (default `false`) and `--save_transcript_secondary_alignment` (default `true`): keep minimap2's secondary alignments. oarfish uses the transcriptome secondary alignments to resolve multi-mapping reads.
- `--kmer_size` (default `14`): the minimap2 k-mer size for the cDNA alignments.

### Parameter files

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](#resource-requests), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run IntGenomicsLab/lrsinglecell -r <version> -profile <docker/singularity/apptainer/institute> -params-file params.yaml
```

with

```yaml title="params.yaml"
input: "./samplesheet.csv"
outdir: "./results/"
genome_fasta: "/path/to/genome.fa"
gtf: "/path/to/file.gtf"
quantifier: "isoquant,oarfish"
barcode_format: "10X_3v3"
```

`quantifier` takes one of `isoquant`, `oarfish` or `lrkallisto`, or a comma-separated combination without spaces such as `isoquant,oarfish,lrkallisto`.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull IntGenomicsLab/lrsinglecell
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [IntGenomicsLab/lrsinglecell releases page](https://github.com/IntGenomicsLab/lrsinglecell/releases) and find the latest pipeline version. Release tags start with a `v` (e.g. `v1.1.0`). Then specify this when running the pipeline with `-r` (one hyphen), e.g. `-r v1.1.0`. You can switch to another version by changing the tag after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports, and in `pipeline_info/`.

To further assist in reproducibility, you can share and reuse [parameter files](#parameter-files) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such a parameter file (for example as supplementary material for a publication), make sure it does not include cluster-specific paths to files, nor institutional profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker, Singularity or Apptainer containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also loads the institutional profiles of [nf-core/configs](https://github.com/nf-core/configs) at run time, so if your cluster has a profile there you can add it, e.g. `-profile singularity,<institute>`. See the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation) for the available profiles. Offline, or for a cluster without one, pass your own config with [`-c`](#-c).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

Test profiles. Each runs the pipeline end to end on small public data from [IntGenomicsLab/test-datasets](https://github.com/IntGenomicsLab/test-datasets/tree/main/lrsinglecell) and needs no other parameters than `--outdir`. All except `test_full` use a chr21 reference and are capped at 4 CPUs and 15 GB.

- `test`
  - cDNA only: subsamples of two 10X 3' v3 libraries (ERR9958133, ERR9958134; `--barcode_format 10X_3v3`), flexiplex demultiplexing, `--quantifier isoquant`.
- `test_dna`
  - DNA only: the DNA library of a 10X Multiome subsample (`--barcode_format 10X_multiome`), flexiplex demultiplexing and DNA duplicate marking.
- `test_mixed`
  - cDNA + DNA: both libraries of the same 10X Multiome subsample (`--barcode_format 10X_multiome`), `--quantifier isoquant`. This is the profile used in the [tutorial](tutorial.md).
- `test_lrkallisto`
  - The cDNA data of `test`, quantified with `--quantifier lrkallisto` (with `--lrkallisto_all_min_reads 0`, so the tiny all-droplet matrix is not filtered).
- `test_oarfish`
  - The cDNA data of `test`, quantified with `--quantifier isoquant,oarfish` and no `--transcript_fasta`, so the oarfish arm aligns to the transcriptome derived from the genome and GTF with `gffread`.
- `test_full`
  - A full-size cDNA run: the complete ERR9958133 and ERR9958134 FASTQs (downloaded from the ENA) against GRCh38 with the GENCODE v45 genome, transcriptome and annotation, `--barcode_format 10X_3v3`, `--split_amount 500000` and all three quantifiers. It has no resource cap and needs a cluster.

Software profiles:

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above.
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). Use it for resource requests, executor and scheduler settings, container options and module arguments, not for pipeline parameters. See the [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time, set by its label in [`conf/base.config`](../conf/base.config) (`process_single`, `process_low`, `process_medium`, `process_high`, ...).

If a task fails with an exit code that usually means it was killed for running out of memory or time (the codes are listed in the `errorStrategy` in `conf/base.config`), it is retried automatically. The pipeline allows two retries (`maxRetries = 2`), so a task gets three attempts in total, and its resource requests are multiplied by the attempt number (`task.attempt`), so the second and third attempts ask for two and three times the resources of the first. Any other failure, or a task that still fails on its last attempt, stops the pipeline once the running tasks have finished.

To change the resource requests of a process, or to cap them at what your machine or cluster offers (`process.resourceLimits`), write a config file and pass it with `-c`. For example:

```groovy title="custom.config"
process {
    resourceLimits = [ cpus: 32, memory: 256.GB, time: 72.h ]

    withName: '.*:ISOQUANT' {
        cpus   = 30
        memory = 85.GB
    }
}
```

See the Nextflow documentation on [process selectors](https://www.nextflow.io/docs/latest/config.html#process-selectors) and [`resourceLimits`](https://www.nextflow.io/docs/latest/reference/process.html#resourcelimits) for the syntax. The execution trace in `pipeline_info/` lists the peak memory (`peak_rss`) and run time of every task, which is the best guide to what needs changing.

### Custom containers

Every process runs in its own container (from [BioContainers](https://biocontainers.pro/), built from [Bioconda](https://bioconda.github.io/) packages) or conda environment, pinned in the module. To use a different container or tool version for a process, override its `container` (or `conda`) directive in a config passed with `-c`:

```groovy title="custom.config"
process {
    withName: '.*:OARFISH' {
        container = 'quay.io/biocontainers/oarfish:<version>'
    }
}
```

Results obtained with a different tool version are no longer those of the pipeline release, so record this when reporting them.

### Custom tool arguments

Most tools receive their arguments from `ext.args` in [`conf/modules.config`](../conf/modules.config). To add or change arguments of a tool, override `ext.args` for its process in a config passed with `-c`, as in the [lr-kallisto strandedness example](#lr-kallisto) above. Copy the existing `ext.args` of the process from `conf/modules.config` first: your value replaces it, it is not appended to it.

### MultiQC report

`--multiqc_config` takes a [MultiQC config file](https://docs.seqera.io/multiqc/getting_started/config) that is applied to both reports (`raw_qc/` and `final_qc/`) after the pipeline's own: MultiQC is given [`assets/multiqc_config.yml`](../assets/multiqc_config.yml), then the report comment with the pipeline version, then your file, each with its own `--config`. A later config wins, so your file only needs the settings you want to change (e.g. `report_header_info`, a `custom_logo` or a different `report_section_order`); the rest of the pipeline's config stays in effect. `--multiqc_logo` adds a logo, and `--multiqc_methods_description` replaces the template of the methods description ([`assets/methods_description_template.yml`](../assets/methods_description_template.yml)).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Troubleshooting

If you experience any issues, please [open an issue on our GitHub repository](https://github.com/IntGenomicsLab/lrsinglecell/issues/new/choose). Some resolutions for common issues are noted below.

### Large runs and `--split_amount`

`--split_amount` (default `0`, no splitting) splits each sample's FASTQ into chunks of that many **reads** (`seqkit split2 --by-size`), so that trimming (chopper) and flexiplex barcode discovery and assignment run on the chunks in parallel. The chunks are merged again before alignment, and the per-chunk flexiplex barcode counts are summed (`flexiplex/<sample>_barcode_counts.txt`). For a PromethION flow cell, `500000` is a good starting value. Do not set it much lower: every chunk is a separate job, so a small value creates a very large number of tasks.

```yml title="params.yml"
split_amount: 500000
```

### Resources for large datasets

The defaults are set for GridION and average-sized PromethION runs. If a task still fails after its automatic retries (see [Resource requests](#resource-requests)), raise its requests in a custom config passed with `-c`. These are overrides we have used on large PromethION datasets; they may not fit every dataset, but they show which processes are the heaviest:

```groovy title="custom.config"
process {
    withName: '.*:.*FASTQC.*' {
        cpus = 20
    }
    withName: '.*:SAMTOOLS_SORT' {
        cpus = 20
    }
    withName: '.*:MINIMAP2_ALIGN' {
        cpus = 20
    }
    withName: '.*:ISOQUANT' {
        cpus   = 30
        memory = 85.GB
    }
}
```

`pipeline_info/execution_trace_<timestamp>.txt` shows the peak memory and run time of each task, which tells you which process to change.

We have also seen recurrent node failures on Slurm clusters that are related to how Nextflow submits jobs rather than to this pipeline. Two things that help: give the failing process more memory with a custom config, or run the whole pipeline inside one interactive session or job with enough CPUs and memory (slower, since there is less parallelisation).

### Transcriptome deduplication and `--fasta_delimiter`

UMI deduplication is split by chromosome for speed. For the transcriptome alignment the pipeline first has to group the transcripts by chromosome, which it does by reading the transcript id from the transcriptome FASTA headers and looking it up in the GTF. Transcript ids are often embedded in a longer, delimited FASTA header; the delimiters used by GENCODE, NCBI and Ensembl references are detected automatically. If your reference uses a different delimiter, set it with `--fasta_delimiter`. This only matters when you give `--transcript_fasta`; the transcriptome derived with `gffread` has plain transcript ids.

### BLAZE

With `--demux_tool_cdna blaze`, barcodes are called by BLAZE and then cut out of the reads and corrected by the pipeline's own scripts, `bin/pre_extract_barcodes.py` (process `PREEXTRACT_FASTQ`) and `bin/correct_barcodes.py` (process `CORRECT_BARCODES`). With `--split_amount` > 0 these two run on chunks of that many reads. Without splitting, `CORRECT_BARCODES` can take hours rather than minutes on a large sample; if you run without `--split_amount`, consider raising its time limit. The BLAZE-only tagging step (`TAG_BARCODES`) holds every read's barcode and UMI in memory and already requests 72 GB and 20 h per attempt; only raise those if its trace shows it running out:

```groovy title="custom.config"
process {
    withName: '.*:BLAZE' {
        cpus = 30
    }
    withName: '.*:CORRECT_BARCODES' {
        time = 15.h
    }
}
```

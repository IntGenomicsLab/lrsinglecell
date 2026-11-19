//
// Subworkflow with functionality specific to the IntGenomicsLab/lrsinglecell pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map{
            meta, fastq, cell_count_val, type ->
                return [ [meta.id, type], meta + [ single_end:true, cell_count: cell_count_val, type: type ], [ fastq ] ]
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()

    def samplesheet_rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
    def has_cdna_input = samplesheet_rows.any { row -> row[3].toString().equalsIgnoreCase('cdna') }
    def has_dna_input  = samplesheet_rows.any { row -> row[3].toString().equalsIgnoreCase('dna') }

    //
    // Demultiplexing tools. The schema enums already reject unknown values;
    // these checks also cover runs with --validate_params false.
    //
    if (!(params.demux_tool_cdna in ['flexiplex', 'blaze'])) {
        error("Unknown --demux_tool_cdna '${params.demux_tool_cdna}'. Please use one of: flexiplex, blaze")
    }
    if (params.demux_tool_dna != 'flexiplex') {
        error("Unknown or unsupported --demux_tool_dna '${params.demux_tool_dna}'. DNA reads can only be demultiplexed with flexiplex")
    }

    //
    // Deduplication
    //
    // The schema enum already rejects anything but umitools; these checks also
    // cover runs with --validate_params false.
    //
    if (params.dedup_tool == 'picard') {
        error("--dedup_tool picard was removed in v1.2.0; use --dedup_tool umitools")
    }
    if (params.dedup_tool != 'umitools') {
        error("Unknown --dedup_tool '${params.dedup_tool}'. Only umitools is supported")
    }

    //
    // DNA input
    //
    if (has_dna_input) {
        if (!params.genome_fasta) {
            error("Input contains DNA reads but --genome_fasta was not provided. DNA reads are aligned to the genome, so a genome fasta is required")
        }

        // Without a preset or a custom pattern flexiplex would search the DNA
        // reads with its RNA default pattern and no whitelist.
        if (params.barcode_format != '10X_multiome' && !params.custom_flexiplex_barcode_dna) {
            error("Input contains DNA reads but no DNA barcode layout is known. Please set --barcode_format 10X_multiome, or give the flexiplex pattern with --custom_flexiplex_barcode_dna (together with --dna_whitelist)")
        }

        if (params.custom_flexiplex_barcode_dna && !params.dna_whitelist) {
            error("--custom_flexiplex_barcode_dna was given without --dna_whitelist. A custom DNA barcode pattern needs a matching barcode whitelist")
        }
    }

    //
    // cDNA input
    //

    // Only require a quantifier when cDNA input is present.
    if (has_cdna_input && !params.quantifier) {
        error("Input contains cDNA reads but --quantifier was not provided. Please set --quantifier to one or more of: isoquant,oarfish,lrkallisto")
    }

    def quantifiers = params.quantifier ? params.quantifier.split(',') as List : []

    if (quantifiers.contains('isoquant') && !params.genome_fasta) {
        error("In order to quantify with isoquant, a genome fasta must be provided")
    }

    if (quantifiers.contains('oarfish')) {
        // Without --transcript_fasta the transcriptome is extracted from the
        // genome and the gtf with gffread (see PREPARE_REFERENCE_FILES)
        if (!params.transcript_fasta && !params.genome_fasta) {
            error("In order to quantify with oarfish, either a transcript fasta (--transcript_fasta) or a genome fasta (--genome_fasta) to extract the transcripts from must be provided")
        }

        // GROUP_TRANSCRIPTS reads the GTF to split the transcriptome by chromosome
        if (!params.gtf) {
            error("In order to quantify with oarfish, a gtf must be provided (--gtf)")
        }
    }

    if (quantifiers.contains('lrkallisto')) {
        if (!params.genome_fasta || !params.gtf) {
            error("In order to quantify with lr-kallisto, a genome fasta and a gtf must be provided")
        }

        // lr-kallisto reads the barcode and umi out of the flexiplex read name;
        // the blaze path keeps them in a separate per-read table instead.
        if (params.demux_tool_cdna != 'flexiplex') {
            error("lr-kallisto quantification requires --demux_tool_cdna flexiplex, but '${params.demux_tool_cdna}' was given")
        }

        // A custom flexiplex pattern can change the barcode and umi widths,
        // which lr-kallisto needs in order to locate them positionally.
        if (params.custom_flexiplex_barcode_cdna) {
            error("lr-kallisto quantification cannot infer the barcode and umi widths from --custom_flexiplex_barcode_cdna. Please use one of the --barcode_format presets instead")
        }
    }

    // Every cDNA arm needs the annotation (gene tagging, isoquant, transcript grouping).
    // DNA-only runs do not, which is why --gtf is not required by the schema.
    if (has_cdna_input && !params.gtf) {
        error("Input contains cDNA reads but --gtf was not provided. A gtf is required to process cDNA reads; only DNA-only runs can omit it")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    // Rows of one sample and type are merged and only the first cell_count would be
    // kept, so they must agree. DNA rows do not use it, but must still be consistent.
    def cell_counts = metas.collect { meta -> meta.cell_count }.unique()
    if (cell_counts.size() != 1) {
        error("Please check input samplesheet -> Rows of sample '${metas[0].id}' (type ${metas[0].type}) have different cell_count values (${cell_counts.join(', ')}). Rows of the same sample and type must use the same cell_count")
    }
    if (cell_counts[0] == null || (cell_counts[0] as Integer) <= 0) {
        error("Please check input samplesheet -> cell_count of sample '${metas[0].id}' (type ${metas[0].type}) must be greater than 0, but is ${cell_counts[0]}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//

//
// In-text citation and bibliography entry of every tool the methods text can name.
// The entries follow CITATIONS.md. NanoPlot and NanoComp share the NanoPack paper,
// and toolBibliographyText lists a shared entry once.
//
def doiLink(id) {
    return "doi: <a href=\"https://doi.org/${id}\">${id}</a>"
}

def toolReferences() {
    return [
        seqkit    : [
            text: 'SeqKit (Shen et al. 2016)',
            bib : "Shen W, Le S, Li Y, Hu F. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS One. 2016 Oct 5;11(10):e0163962. ${doiLink('10.1371/journal.pone.0163962')}"
        ],
        chopper   : [
            text: 'chopper (De Coster and Rademakers 2023)',
            bib : "De Coster W, Rademakers R. NanoPack2: population-scale evaluation of long-read sequencing data. Bioinformatics. 2023 May;39(5):btad311. ${doiLink('10.1093/bioinformatics/btad311')}"
        ],
        flexiplex : [
            text: 'flexiplex (Cheng et al. 2024)',
            bib : "Cheng O, Ling MH, Wang C, Wu S, Ritchie ME, Göke J, Amin N, Davidson NM. Flexiplex: a versatile demultiplexer and search tool for omics data. Bioinformatics. 2024 Mar 4;40(3):btae102. ${doiLink('10.1093/bioinformatics/btae102')}"
        ],
        blaze     : [
            text: 'BLAZE (You et al. 2023)',
            bib : "You Y, Prawer YDJ, De Paoli-Iseppi R, Hunt CPJ, Parish CL, Shim H, Clark MB. Identification of cell barcodes from long-read single-cell RNA-seq with BLAZE. Genome Biol. 2023 Apr 6;24(1):66. ${doiLink('10.1186/s13059-023-02907-y')}"
        ],
        minimap2  : [
            text: 'minimap2 (Li 2018)',
            bib : "Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018 Sep 15;34(18):3094-3100. ${doiLink('10.1093/bioinformatics/bty191')}"
        ],
        samtools  : [
            text: 'SAMtools (Li et al. 2009)',
            bib : "Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. ${doiLink('10.1093/bioinformatics/btp352')}"
        ],
        pysam     : [
            text: 'pysam',
            bib : 'pysam: a Python module for reading and manipulating SAM/BAM/CRAM files. URL: <a href="https://github.com/pysam-developers/pysam">https://github.com/pysam-developers/pysam</a>'
        ],
        bamtools  : [
            text: 'BamTools (Barnett et al. 2011)',
            bib : "Barnett DW, Garrison EK, Quinlan AR, Strömberg MP, Marth GT. BamTools: a C++ API and toolkit for analyzing and managing BAM files. Bioinformatics. 2011 Jun 15;27(12):1691-2. ${doiLink('10.1093/bioinformatics/btr174')}"
        ],
        umitools  : [
            text: 'UMI-tools (Smith et al. 2017)',
            bib : "Smith T, Heger A, Sudbery I. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Res. 2017 Mar;27(3):491-499. ${doiLink('10.1101/gr.209601.116')}"
        ],
        isoquant  : [
            text: 'IsoQuant (Prjibelski et al. 2023)',
            bib : "Prjibelski AD, Mikheenko A, Joglekar A, Smetanin A, Jarroux J, Lapidus AL, Tilgner HU. Accurate isoform discovery with IsoQuant using long reads. Nat Biotechnol. 2023 Jul;41(7):915-918. ${doiLink('10.1038/s41587-022-01565-y')}"
        ],
        oarfish   : [
            text: 'oarfish (Jousheghani et al. 2025)',
            bib : "Jousheghani ZZ, Singh NP, Patro R. Oarfish: Enhanced probabilistic modeling leads to improved accuracy in long read transcriptome quantification. Bioinformatics. 2025 Jul 1;41(Supplement_1):i304-i313. ${doiLink('10.1093/bioinformatics/btaf240')}"
        ],
        lrkallisto: [
            text: 'lr-kallisto (Loving et al. 2025)',
            bib : "Loving RK, Sullivan DK, Reese F, Rebboah E, Sakr J, Rezaie N, Liang HY, Filimban G, Kawauchi S, Booeshaghi AS, Melsted P, Oakes C, Trout D, Williams BA, MacGregor GR, Wold BJ, Mortazavi A, Pachter L. Long-read sequencing transcriptome quantification with lr-kallisto. PLoS Comput Biol. 2025 Dec 1;21(12):e1013692. ${doiLink('10.1371/journal.pcbi.1013692')}"
        ],
        kallisto  : [
            text: 'kallisto (Bray et al. 2016)',
            bib : "Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol. 2016 May;34(5):525-7. ${doiLink('10.1038/nbt.3519')}"
        ],
        bustools  : [
            text: 'bustools (Melsted et al. 2021)',
            bib : "Melsted P, Booeshaghi AS, Liu L, Gao F, Lu L, Min KHJ, da Veiga Beltrame E, Hjörleifsson KE, Gehring J, Pachter L. Modular, efficient and constant-memory single-cell RNA-seq preprocessing. Nat Biotechnol. 2021 Jul;39(7):813-818. ${doiLink('10.1038/s41587-021-00870-2')}"
        ],
        gffread   : [
            text: 'gffread (Pertea and Pertea 2020)',
            bib : "Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. F1000Res. 2020 Apr 28;9:ISCB Comm J-304. ${doiLink('10.12688/f1000research.23297.2')}"
        ],
        nanoplot  : [
            text: 'NanoPlot (De Coster et al. 2018)',
            bib : "De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018 Aug 1;34(15):2666-2669. ${doiLink('10.1093/bioinformatics/bty149')}"
        ],
        nanocomp  : [
            text: 'NanoComp (De Coster et al. 2018)',
            bib : "De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018 Aug 1;34(15):2666-2669. ${doiLink('10.1093/bioinformatics/bty149')}"
        ],
        toulligqc : [
            text: 'ToulligQC (Dias et al.)',
            bib : 'Dias K, Laffay B, Ferrato-Berberian L, Lemoine S, Hamraoui A, Thomas-Chollier M, Le Crom S, Jourdren L. ToulligQC. URL: <a href="https://github.com/GenomiqueENS/toulligQC">https://github.com/GenomiqueENS/toulligQC</a>'
        ],
        fastqc    : [
            text: 'FastQC (Andrews 2010)',
            bib : 'Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. URL: <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a>'
        ],
        rseqc     : [
            text: 'RSeQC (Wang et al. 2012)',
            bib : "Wang L, Wang S, Li W. RSeQC: quality control of RNA-seq experiments. Bioinformatics. 2012 Aug 15;28(16):2184-5. ${doiLink('10.1093/bioinformatics/bts356')}"
        ],
        ucsc      : [
            text: 'the UCSC tools gtfToGenePred and genePredToBed (Kent et al. 2010)',
            bib : "Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics. 2010 Sep 1;26(17):2204-7. ${doiLink('10.1093/bioinformatics/btq351')}"
        ],
        seurat    : [
            text: 'Seurat (Hao et al. 2021)',
            bib : "Hao Y, Hao S, Andersen-Nissen E, Mauck WM 3rd, Zheng S, Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, Hoffman P, Stoeckius M, Papalexi E, Mimitou EP, Jain J, Srivastava A, Stuart T, Fleming LM, Yeung B, Rogers AJ, McElrath JM, Blish CA, Gottardo R, Smibert P, Satija R. Integrated analysis of multimodal single-cell data. Cell. 2021 Jun 24;184(13):3573-3587.e29. ${doiLink('10.1016/j.cell.2021.04.048')}"
        ],
        multiqc   : [
            text: 'MultiQC (Ewels et al. 2016)',
            bib : "Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. ${doiLink('10.1093/bioinformatics/btw354')}"
        ],
    ]
}

//
// The steps this run performs, in pipeline order, as [ sentence, [ tool keys ] ].
// Each step is included only when the parameters, and for the cDNA and DNA arms the
// sample types in the samplesheet (sample_types: e.g. [ 'cdna', 'dna' ]), mean it
// runs, so the text describes this run.
//
def methodsSteps(sample_types) {
    def refs = toolReferences()

    def has_cdna = sample_types.contains('cdna')
    def has_dna  = sample_types.contains('dna')

    def quants            = params.quantifier ? params.quantifier.tokenize(',')*.trim() : []
    def genome_arm        = has_cdna && (quants.contains('isoquant') || quants.contains('lrkallisto'))
    def transcriptome_arm = has_cdna && quants.contains('oarfish')
    def run_qc            = !params.skip_qc

    def steps = []

    // Pre-processing
    // FILTER_DUPLEX_PARENTS is the pipeline's own awk filter, so there is no tool to cite
    if (params.drop_duplex_parents.toString().toBoolean()) {
        steps << [ "The simplex parent reads (dx:i:-1) of dorado duplex reads were removed, so that each duplex molecule was used once, as its consensus read.", [] ]
    }
    // SEQKIT_SPLIT2 only runs as part of trimming (workflows/lrsinglecell.nf)
    if (!params.skip_trimming && (params.split_amount as Integer) > 0) {
        steps << [ "Reads were split into chunks for parallel processing with ${citeTools(refs, ['seqkit'])}.", [ 'seqkit' ] ]
    }
    if (!params.skip_trimming) {
        steps << [ "Reads were quality-filtered and trimmed with ${citeTools(refs, ['chopper'])}.", [ 'chopper' ] ]
    }

    // Demultiplexing
    if (has_cdna) {
        def demux = params.demux_tool_cdna == 'blaze' ? 'blaze' : 'flexiplex'
        steps << [ "Cell barcodes and UMIs of the cDNA reads were identified with ${citeTools(refs, [demux])}.", [ demux ] ]
    }
    if (has_dna) {
        steps << [ "Cell barcodes of the DNA reads were identified with ${citeTools(refs, ['flexiplex'])}.", [ 'flexiplex' ] ]
    }

    // Alignment
    def cdna_refs = []
    if (genome_arm)        { cdna_refs << 'the genome' }
    if (transcriptome_arm) { cdna_refs << 'the transcriptome' }
    def aligned = []
    if (cdna_refs) { aligned << "cDNA reads were aligned to ${cdna_refs.join(' and ')}" }
    if (has_dna)   { aligned << (aligned ? 'DNA reads to the genome' : 'DNA reads were aligned to the genome') }
    if (aligned) {
        steps << [
            "${aligned.join(' and ')} with ${citeTools(refs, ['minimap2'])}. Alignments were sorted, indexed, filtered and summarised with ${citeTools(refs, ['samtools'])} and ${citeTools(refs, ['pysam'])}.",
            [ 'minimap2', 'samtools', 'pysam' ]
        ]
    }

    // Deduplication
    def genome_dedup = genome_arm && !params.skip_cdna_dedup
    if (params.dedup_tool == 'umitools' && (genome_dedup || transcriptome_arm)) {
        def keys = [ 'umitools' ]
        def detail = ''
        if (genome_dedup) {
            keys << 'bamtools'
            detail = params.dedup_per_gene
                ? "; genome alignments were split per chromosome with ${citeTools(refs, ['bamtools'])}, assigned to genes by gene-body overlap with ${citeTools(refs, ['pysam'])} and grouped by gene"
                : "; genome alignments were split per chromosome with ${citeTools(refs, ['bamtools'])}"
        }
        steps << [ "cDNA UMI duplicates were removed with ${citeTools(refs, ['umitools'])}${detail}.", keys ]
    }
    if (has_dna && !params.skip_dedup) {
        steps << [ "Duplicate DNA reads were marked per cell barcode with the pipeline's mark_dna_duplicates.py, built on ${citeTools(refs, ['pysam'])}.", [ 'pysam' ] ]
    }

    // Quantification
    if (has_cdna && quants.contains('isoquant')) {
        steps << [ "Genes and transcripts were quantified per cell from the genome alignments with ${citeTools(refs, ['isoquant'])}.", [ 'isoquant' ] ]
    }
    if (transcriptome_arm) {
        // Without --transcript_fasta the transcriptome is extracted from the genome
        // and the gtf with gffread (PREPARE_REFERENCE_FILES)
        def derived = params.transcript_fasta
            ? ''
            : ", using transcript sequences extracted from the genome and the GTF with ${citeTools(refs, ['gffread'])}"
        steps << [
            "Transcripts were quantified per cell from the transcriptome alignments with ${citeTools(refs, ['oarfish'])}${derived}.",
            params.transcript_fasta ? [ 'oarfish' ] : [ 'oarfish', 'gffread' ]
        ]
    }
    if (has_cdna && quants.contains('lrkallisto')) {
        steps << [
            "Genes and transcripts were quantified per cell with ${citeTools(refs, ['lrkallisto'])}, using ${citeToolsAnd(refs, ['kallisto', 'bustools'])}, on transcript sequences extracted from the GTF with ${citeTools(refs, ['gffread'])}.",
            [ 'lrkallisto', 'kallisto', 'bustools', 'gffread' ]
        ]
    }

    // Quality control
    def qc = []
    if (run_qc) {
        if (!params.skip_nanoplot)                                    { qc << 'nanoplot' }
        if (!params.skip_toulligqc)                                   { qc << 'toulligqc' }
        if (!params.skip_fastqc)                                      { qc << 'fastqc' }
        if (!params.skip_fastq_nanocomp || !params.skip_bam_nanocomp) { qc << 'nanocomp' }
    }
    if (qc) {
        steps << [ "Read quality was assessed with ${citeToolsAnd(refs, qc)}.", qc ]
    }
    if (run_qc && !params.skip_rseqc && params.gtf && genome_arm) {
        steps << [
            "The distribution of cDNA reads over genomic features was computed with ${citeTools(refs, ['rseqc'])}, on a BED file made from the GTF with ${citeTools(refs, ['ucsc'])}.",
            [ 'rseqc', 'ucsc' ]
        ]
    }
    if (run_qc && !params.skip_seurat && has_cdna && quants) {
        steps << [ "Per-cell summary statistics of the unfiltered count matrices were computed with ${citeTools(refs, ['seurat'])}.", [ 'seurat' ] ]
    }
    steps << [ "Results were summarised in a report with ${citeTools(refs, ['multiqc'])}.", [ 'multiqc' ] ]

    return steps
}

// "A (ref), B (ref)": the in-text citations of the given tools
def citeTools(refs, keys) {
    return keys.collect { key -> refs[key].text }.join(', ')
}

// "A (ref), B (ref) and C (ref)"
def citeToolsAnd(refs, keys) {
    return keys.size() > 1 ? citeTools(refs, keys[0..-2]) + ' and ' + refs[keys[-1]].text : citeTools(refs, keys)
}

def toolCitationText(sample_types) {
    return methodsSteps(sample_types).collect { step -> step[0] }.join(' ')
}

def toolBibliographyText(sample_types) {
    def refs = toolReferences()
    return methodsSteps(sample_types)
        .collectMany { step -> step[1] }
        .collect { key -> refs[key].bib.toString() }
        .unique()
        .collect { bib -> "<li>${bib}</li>" }
        .join('\n    ')
}

//
// The MultiQC report_comment, linking to the exact pipeline version that ran: the
// release tag for a release, the dev branch for a development version. The release
// tags of this repository carry a "v" prefix (v1.1.0), so this is built at run time
// rather than left to `nf-core pipelines bump-version`, which writes the version
// without it. It is passed to MultiQC as an extra config after assets/multiqc_config.yml.
//
def multiqcReportComment() {
    def manifest = workflow.manifest
    def version  = manifest.version ?: ''
    def name     = manifest.name ?: 'IntGenomicsLab/lrsinglecell'
    def home     = manifest.homePage ?: "https://github.com/${name}"
    def ref      = (!version || version.contains('dev'))
        ? 'tree/dev'
        : "releases/tag/${version.startsWith('v') ? version : 'v' + version}"
    def label    = version ? " v${version.replaceFirst(/^v/, '')}" : ''
    return "report_comment: >\n" +
        "  This report has been generated by the <a href=\"${home}/${ref}\" target=\"_blank\">${name}${label}</a> analysis pipeline.\n"
}

def methodsDescriptionText(mqc_methods_yaml, sample_types) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of the version of the pipeline used. </li>"

    // Tool references, for the steps this run performs
    meta["tool_citations"] = toolCitationText(sample_types)
    meta["tool_bibliography"] = toolBibliographyText(sample_types)

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

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
    if (params.dedup_tool == 'picard') {
        log.warn("--dedup_tool picard is deprecated and will be removed in v1.1.0. It applies to cDNA only (DNA is always deduplicated by mark_dna_duplicates.py), groups on the cell barcode and ignores the UMI. Please use --dedup_tool umitools instead")
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
        if (!params.transcript_fasta) {
            error("In order to quantify with oarfish, a transcript fasta must be provided")
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
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
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
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

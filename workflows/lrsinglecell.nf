/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { CHOPPER                           } from "../modules/local/chopper"
include { FILTER_DUPLEX_PARENTS             } from "../modules/local/filter_duplex_parents"
include { READ_COUNTS                       } from "../modules/local/read_counts"
include { UCSC_GENEPREDTOBED                } from "../modules/local/ucsc_genepredtobed"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { PREPARE_REFERENCE_FILES                                     } from "../subworkflows/local/prepare_reference_files"
include { DEMULTIPLEX_FLEXIPLEX as DEMULTIPLEX_FLEXIPLEX_CDNA         } from "../subworkflows/local/demultiplex_flexiplex"
include { DEMULTIPLEX_FLEXIPLEX as DEMULTIPLEX_FLEXIPLEX_DNA          } from "../subworkflows/local/demultiplex_flexiplex"
include { DEMULTIPLEX_BLAZE                                           } from "../subworkflows/local/demultiplex_blaze"
include { PROCESS_LONGREAD_SCRNA as PROCESS_LONGREAD_SCRNA_GENOME     } from "../subworkflows/local/process_longread_scrna"
include { PROCESS_LONGREAD_SCRNA as PROCESS_LONGREAD_SCRNA_TRANSCRIPT } from "../subworkflows/local/process_longread_scrna"
include { QUANTIFY_SCRNA_LRKALLISTO                                   } from "../subworkflows/local/quantify_scrna_lrkallisto"
include { ALIGN_DEDUPLICATE_DNA                                       } from "../subworkflows/local/align_deduplicate_dna"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Installed directly from nf-core/modules
//
include { NANOCOMP as NANOCOMP_FASTQ_CDNA               } from "../modules/nf-core/nanocomp/main"
include { NANOCOMP as NANOCOMP_FASTQ_DNA                } from "../modules/nf-core/nanocomp/main"
include { SEQKIT_SPLIT2                                 } from "../modules/nf-core/seqkit/split2"
include { MULTIQC as MULTIQC_RAWQC                      } from "../modules/nf-core/multiqc/main"
include { MULTIQC as MULTIQC_FINALQC                    } from "../modules/nf-core/multiqc/main"
include { CAT_FASTQ                                     } from "../modules/nf-core/cat/fastq/main"
include { CAT_CAT                                       } from "../modules/nf-core/cat/cat/main"
include { UCSC_GTFTOGENEPRED                            } from "../modules/nf-core/ucsc/gtftogenepred/main"
include { paramsSummaryMap                              } from "plugin/nf-schema"

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/subworkflows
 */

include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_PRE_TRIM          } from "../subworkflows/nf-core/qcfastq_nanoplot_fastqc"
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_TRIM         } from "../subworkflows/nf-core/qcfastq_nanoplot_fastqc"
include { QCFASTQ_NANOPLOT_FASTQC as FASTQC_NANOPLOT_POST_EXTRACT      } from "../subworkflows/nf-core/qcfastq_nanoplot_fastqc"
include { paramsSummaryMultiqc                                         } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { softwareVersionsToYAML                                       } from "../subworkflows/nf-core/utils_nfcore_pipeline"
include { methodsDescriptionText                                       } from "../subworkflows/local/utils_nfcore_lrsinglecell_pipeline"
include { multiqcReportComment                                         } from "../subworkflows/local/utils_nfcore_lrsinglecell_pipeline"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LRSINGLECELL {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PARAMETER PRESETS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Whitelist
    def cdna_whitelist = null
    if (params.cdna_whitelist) {
        cdna_whitelist = file(params.cdna_whitelist)
    }
    else if (params.barcode_format.equals("10X_3v3")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/3M-february-2018.txt.gz")
    }
    else if (params.barcode_format.equals("10X_5v2")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/737K-august-2016.txt.gz")
    }
    else if (params.barcode_format.equals("10X_3v4")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/3M-3pgex-may-2023_TRU.txt.gz")
    }
    else if (params.barcode_format.equals("10X_5v3")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/3M-5pgex-jan-2023.txt.gz")
    }
    else if (params.barcode_format.equals("10X_multiome")) {
        cdna_whitelist = file("$baseDir/assets/whitelist/cellranger_arc_rna.737K-arc-v1.txt.gz")
    }
    else {
        cdna_whitelist = []
    }

    def dna_whitelist = null
    if (params.dna_whitelist) {
        dna_whitelist = file(params.dna_whitelist)
    }
    else if (params.barcode_format.equals("10X_multiome")) {
        dna_whitelist = file("$baseDir/assets/whitelist/cellranger_arc_atac.737K-arc-v1.txt.gz")
    }
    else {
        dna_whitelist = []
    }

    // Quantifiers
    // Associate the quantifiers with the kind of alignment needed
    def GENOME_QUANT_OPTS     = [ 'isoquant' ]
    def TRANSCRIPT_QUANT_OPTS = [ 'oarfish' ]
    // lr-kallisto pseudoaligns straight from the demultiplexed fastq, so it
    // needs no alignment and no separate deduplication step
    def FASTQ_QUANT_OPTS      = [ 'lrkallisto' ]

    def quantifier_list   = params.quantifier ? params.quantifier.split(',') as List : []
    def genome_quants     = quantifier_list.findAll { q -> q in GENOME_QUANT_OPTS }
    def transcript_quants = quantifier_list.findAll { q -> q in TRANSCRIPT_QUANT_OPTS }
    def fastq_quants      = quantifier_list.findAll { q -> q in FASTQ_QUANT_OPTS }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONFIG FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def ch_multiqc_config                     = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config              = multiqc_config ? channel.fromPath( multiqc_config, checkIfExists: true ) : channel.empty()
    def ch_multiqc_logo                       = multiqc_logo   ? channel.fromPath( multiqc_logo, checkIfExists: true ) : channel.empty()
    def ch_multiqc_custom_methods_description = multiqc_methods_description ? channel.fromPath(multiqc_methods_description, checkIfExists: true) : channel.fromPath("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    def ch_versions = channel.empty()
    def ch_multiqc_report = channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    ch_samplesheet
        .branch {
            meta, fastq ->
                single: fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastqs }


    //
    // MODULE: Combine fastqs from the same sample
    //

    CAT_FASTQ ( ch_fastqs.multiple )
        .reads
        .mix ( ch_fastqs.single )
        .set { ch_cat_fastq }

    //
    // MODULE: Drop the simplex parents (dx:i:-1) of dorado duplex reads, so each
    // duplex molecule enters once, as its consensus read. FASTQs without dx tags
    // pass through unchanged. Runs before QC so all counts describe the reads used.
    //

    // `--drop_duplex_parents false` on the command line arrives as the string
    // "false", which is truthy in Groovy, so coerce it rather than testing it directly
    if (params.drop_duplex_parents.toString().toBoolean()) {
        FILTER_DUPLEX_PARENTS ( ch_cat_fastq )
        ch_cat_fastq = FILTER_DUPLEX_PARENTS.out.reads
        ch_versions = ch_versions.mix(FILTER_DUPLEX_PARENTS.out.versions_filter_duplex_parents)

        // A duplex FASTQ without dx tags fails in the module itself; warn about
        // the files that went through unfiltered
        FILTER_DUPLEX_PARENTS.out.mode
            .subscribe { meta, mode ->
                if (mode == 'passthrough_no_dx') {
                    log.warn("--drop_duplex_parents: no dorado 'dx' tags in the reads of ${meta.id} (${meta.type}), so no duplex parents were dropped. " +
                        "This is expected for simplex FASTQs. For duplex data, convert dorado's duplex BAM with 'samtools fastq -T \"*\"' to keep the tags")
                } else if (mode == 'kept_all_no_consensus') {
                    log.warn("--drop_duplex_parents: the reads of ${meta.id} (${meta.type}) have duplex parents (dx:i:-1) but no duplex consensus reads (dx:i:1). " +
                        "Kept all reads, because dropping the parents would lose these molecules")
                }
            }
    }

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot, ToulligQC and FastQC - pre-trim QC
    //

    ch_fastqc_multiqc_pretrim = channel.empty()
    if (!params.skip_qc){

        FASTQC_NANOPLOT_PRE_TRIM ( ch_cat_fastq, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

        ch_fastqc_multiqc_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.fastqc_multiqc.ifEmpty([])
        ch_nanostat_pretrim = FASTQC_NANOPLOT_PRE_TRIM.out.nanoplot_txt.ifEmpty([])
    }

    //
    // MODULE: NanoComp for FastQ files
    //

    ch_nanocomp_fastq_html = channel.empty()
    ch_nanocomp_fastq_txt = channel.empty()
    if (!params.skip_qc && !params.skip_fastq_nanocomp) {

        NANOCOMP_FASTQ_CDNA (
            ch_cat_fastq
                .filter{ meta, fastq -> meta.type == 'cdna' }
                .collect{it[1]}
                .map{
                    [ [ 'id': 'cdna_fastq.', 'type': 'cdna' ] , it ]
                }
        )

        NANOCOMP_FASTQ_DNA (
            ch_cat_fastq
                .filter{ meta, fastq -> meta.type == 'dna' }
                .collect{it[1]}
                .map{
                    [ [ 'id': 'dna_fastq.', 'type': 'dna' ] , it ]
                }
        )

        ch_nanocomp_fastq_html = NANOCOMP_FASTQ_CDNA.out.report_html.mix( NANOCOMP_FASTQ_DNA.out.report_html )
        ch_nanocomp_fastq_txt  = NANOCOMP_FASTQ_CDNA.out.stats_txt.mix( NANOCOMP_FASTQ_DNA.out.stats_txt )

        ch_versions = ch_versions.mix( NANOCOMP_FASTQ_CDNA.out.versions )
        ch_versions = ch_versions.mix( NANOCOMP_FASTQ_DNA.out.versions )

    }

    //
    // SUBWORKFLOW: Prepare reference files
    //

    PREPARE_REFERENCE_FILES (
        params.genome_fasta,
        params.transcript_fasta,
        params.gtf,
        quantifier_list
    )

    genome_fasta = PREPARE_REFERENCE_FILES.out.prepped_genome_fasta
    genome_fai = PREPARE_REFERENCE_FILES.out.genome_fai
    transcript_fasta = PREPARE_REFERENCE_FILES.out.prepped_transcript_fasta
    transcript_fai = PREPARE_REFERENCE_FILES.out.transcript_fai
    gtf = PREPARE_REFERENCE_FILES.out.prepped_gtf

    //
    // MODULE: Generate bed file from input gtf for rseqc
    //

    // come back to this once intron work is finished (likely input will be fine)
    ch_pred = channel.empty()
    ch_rseqc_bed = channel.empty()

    // The bed is only consumed by RSeQC on the cDNA genome arm, which runs for
    // genome and fastq quantifiers. DNA-only runs skip it, and may have no gtf.
    if (!params.skip_qc && !params.skip_rseqc && params.gtf && (genome_quants || fastq_quants)) {
        UCSC_GTFTOGENEPRED( gtf )
        ch_pred = UCSC_GTFTOGENEPRED.out.genepred

        UCSC_GENEPREDTOBED ( ch_pred.map{ _meta, genepred -> [genepred]} )
        ch_rseqc_bed = UCSC_GENEPREDTOBED.out.bed
        ch_versions = ch_versions.mix(UCSC_GENEPREDTOBED.out.versions_genepredtobed)
    }

    //
    // MODULE: Trim and filter reads
    //
    ch_fastqc_multiqc_postrim = channel.empty()
    ch_nanostat_posttrim = channel.empty()
    ch_trimmed_reads_combined = channel.empty()

    if (!params.skip_trimming){

        //
        // MODULE: Chopper
        //

        if ((params.split_amount as Integer) > 0) {
            SEQKIT_SPLIT2( ch_cat_fastq )

            // Temporarily change the meta object so that the id is present on the
            // fastq to prevent duplicated names
            SEQKIT_SPLIT2.out.reads
                .transpose()
                .set { ch_fastqs }

        } else {
            ch_fastqs = ch_cat_fastq
        }

        CHOPPER ( ch_fastqs, params.split_amount )
        ch_trimmed_reads = CHOPPER.out.reads
        ch_versions = ch_versions.mix(CHOPPER.out.versions_chopper)

        // If the fastqs were split, combine them together
        if ((params.split_amount as Integer) > 0){
            CAT_CAT(ch_trimmed_reads.groupTuple())
            ch_trimmed_reads_combined = CAT_CAT.out.file_out
        } else {
            ch_trimmed_reads_combined = ch_trimmed_reads
        }

        //
        // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-trim QC
        //
        if (!params.skip_qc){
            //
            // MODULE: Run qc on the post trimmed reads
            //
            FASTQC_NANOPLOT_POST_TRIM ( ch_trimmed_reads_combined, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

            ch_fastqc_multiqc_postrim = FASTQC_NANOPLOT_POST_TRIM.out.fastqc_multiqc.ifEmpty([])
            ch_nanostat_posttrim = FASTQC_NANOPLOT_POST_TRIM.out.nanoplot_txt.ifEmpty([])
        }
    } else {
        ch_trimmed_reads_combined = ch_cat_fastq
    }


    // Branch channel to dna and cdna
    ch_trimmed_reads_combined = ch_trimmed_reads_combined
        .branch {
            meta, fastq ->
                dna: meta.type == 'dna'
                    return [ meta, fastq ]
                cdna: meta.type == 'cdna'
                    return [ meta, fastq ]
        }

    //
    // SUBWORKFLOW: Demultiplex reads using FLEXIPLEX for DNA
    //

    ch_extracted_fastq_dna = channel.empty()
    ch_corrected_bc_info_dna = channel.empty()
    if (params.demux_tool_dna == "flexiplex") {
        DEMULTIPLEX_FLEXIPLEX_DNA (
            ch_trimmed_reads_combined.dna,
            dna_whitelist
        )

        ch_versions = ch_versions.mix(DEMULTIPLEX_FLEXIPLEX_DNA.out.versions)
        ch_extracted_fastq_dna = DEMULTIPLEX_FLEXIPLEX_DNA.out.flexiplex_fastq
        ch_corrected_bc_info_dna = DEMULTIPLEX_FLEXIPLEX_DNA.out.flexiplex_barcodes
    }
    // Any other --demux_tool_dna is rejected by the schema and validateInputParameters.

    ch_extracted_fastq_cdna = channel.empty()
    ch_corrected_bc_info_cdna = channel.empty()
    // Every barcode flexiplex saw with its read count, before the knee call.
    // Only the flexiplex path produces it; blaze leaves it empty, which is fine
    // because lr-kallisto refuses blaze anyway (see the pipeline validation).
    ch_barcode_counts_cdna = channel.empty()
    if (params.demux_tool_cdna == "flexiplex") {

        //
        // SUBWORKFLOW: Demultiplex reads using FLEXIPLEX for cDNA
        //

        DEMULTIPLEX_FLEXIPLEX_CDNA (
            ch_trimmed_reads_combined.cdna,
            cdna_whitelist
        )

        ch_versions = ch_versions.mix(DEMULTIPLEX_FLEXIPLEX_CDNA.out.versions)

        ch_extracted_fastq_cdna = DEMULTIPLEX_FLEXIPLEX_CDNA.out.flexiplex_fastq
        ch_corrected_bc_info_cdna = DEMULTIPLEX_FLEXIPLEX_CDNA.out.flexiplex_barcodes
        ch_barcode_counts_cdna = DEMULTIPLEX_FLEXIPLEX_CDNA.out.flexiplex_barcode_counts

    } else if (params.demux_tool_cdna == "blaze") {

        //
        // SUBWORKFLOW: Demultiplex reads using BLAZE for cDNA
        //

        DEMULTIPLEX_BLAZE (
            ch_trimmed_reads_combined.cdna,
            cdna_whitelist
        )

        ch_versions = ch_versions.mix(DEMULTIPLEX_BLAZE.out.versions)

        ch_extracted_fastq_cdna = DEMULTIPLEX_BLAZE.out.extracted_fastq
        ch_corrected_bc_info_cdna = DEMULTIPLEX_BLAZE.out.corrected_bc_info
    }

    // Recombine channels for QC modules
    ch_extracted_fastq = ch_extracted_fastq_cdna.mix(ch_extracted_fastq_dna)
    ch_corrected_bc_info = ch_corrected_bc_info_cdna.mix(ch_corrected_bc_info_dna)

    //
    // SUBWORKFLOW: Fastq QC with Nanoplot and FastQC - post-extract QC
    //
    ch_fastqc_multiqc_postextract = channel.empty()
    ch_read_counts = channel.empty()
    if (!params.skip_qc){
        FASTQC_NANOPLOT_POST_EXTRACT ( ch_extracted_fastq, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

        ch_fastqc_multiqc_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.fastqc_multiqc.ifEmpty([])
        ch_nanostat_postextract = FASTQC_NANOPLOT_POST_EXTRACT.out.nanoplot_txt.ifEmpty([])

        //
        // MODULE: Generate read counts
        //

        ch_pretrim_counts = channel.empty()
        ch_posttrim_counts = channel.empty()
        ch_postextract_counts = channel.empty()
        if (!params.skip_fastqc){
            ch_pretrim_counts = ch_fastqc_multiqc_pretrim.collect{it -> it[0]}
            ch_posttrim_counts = ch_fastqc_multiqc_postrim.collect{it -> it[0]}
            ch_postextract_counts = ch_fastqc_multiqc_postextract.collect{it -> it[0]}

        } else if (!params.skip_nanoplot){
            ch_pretrim_counts = ch_nanostat_pretrim.collect{it -> it[1]}
            ch_posttrim_counts = ch_nanostat_posttrim.collect{it -> it[1]}
            ch_postextract_counts = ch_nanostat_postextract.collect{it -> it[1]}

        }

        READ_COUNTS (
            ch_pretrim_counts.ifEmpty([]),
            ch_posttrim_counts.ifEmpty([]),
            ch_postextract_counts.ifEmpty([]),
            ch_corrected_bc_info.collect{it -> it[1]})

        ch_read_counts = READ_COUNTS.out.read_counts
        ch_versions = ch_versions.mix(READ_COUNTS.out.versions_read_counts)
    }

    //
    // SUBWORKFLOW: Align Long Read cDNA data
    //

    ch_multiqc_finalqc_files = channel.empty()

    // Read counts cover every sample, cDNA and DNA, so they are added once here
    // rather than by each quantifier block.
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ch_read_counts.collect().ifEmpty([])
    )

    // The genome alignment is produced whenever there are cDNA reads to align,
    // not only when a genome quantifier asked for one: the bam is a deliverable
    // in its own right (coverage, igv, fusions, variants). lr-kallisto
    // pseudoaligns from the fastq and consumes nothing from here, so with only
    // fastq quantifiers requested this runs the aligner and the deduplication
    // and no quantifier at all -- PROCESS_LONGREAD_SCRNA guards each of those on
    // val_quant_list.contains(...), so an empty genome_quants is safe.
    if (genome_quants || fastq_quants){
        PROCESS_LONGREAD_SCRNA_GENOME(
            genome_fasta,
            genome_fai,
            gtf,
            ch_extracted_fastq_cdna,
            ch_rseqc_bed,
            ch_corrected_bc_info_cdna,
            genome_quants,
            params.dedup_tool,
            true, // Used to indicate the bam is genome aligned
            params.fasta_delimiter,
            params.skip_save_minimap2_index,
            params.skip_qc,
            params.skip_rseqc,
            params.skip_bam_nanocomp,
            params.skip_seurat,
            // Deliberately NOT params.skip_dedup: that one also governs the DNA
            // arm's duplicate marking below, so a single flag could not turn
            // the expensive cDNA umi_tools chain down without silently
            // un-deduplicating the DNA bam too.
            params.skip_cdna_dedup
        )
        ch_versions = ch_versions.mix(PROCESS_LONGREAD_SCRNA_GENOME.out.versions)

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_flagstat.collect{it -> it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_idxstats.collect{it -> it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_rseqc_read_dist.collect{it -> it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.minimap_nanocomp_bam_txt.collect{it -> it[1]}.ifEmpty([])
        )

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.bc_tagged_flagstat.collect{it -> it[1]}.ifEmpty([])
        )

        // With --skip_cdna_dedup, dedup_flagstat is the tagged flagstat above, so
        // gate on the same flag that turns the cDNA dedup off to avoid adding it twice.
        if (!params.skip_cdna_dedup) {
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
                PROCESS_LONGREAD_SCRNA_GENOME.out.dedup_flagstat.collect{it -> it[1]}.ifEmpty([])
            )
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
                PROCESS_LONGREAD_SCRNA_GENOME.out.dedup_idxstats.collect{it -> it[1]}.ifEmpty([])
            )
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
                PROCESS_LONGREAD_SCRNA_GENOME.out.dedup_summary.collect{it -> it[1]}.ifEmpty([])
            )
            // Share of alignments TAG_GENES could not place on a gene, shown in the
            // general statistics as "% reads without gene". Only written with
            // --dedup_per_gene, so empty otherwise.
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
                PROCESS_LONGREAD_SCRNA_GENOME.out.gene_tag_summary.collect{it -> it[1]}.ifEmpty([])
            )
        }

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.gene_qc_stats
                .map{ _meta, stats -> [stats]}
                .collect()
                .ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_GENOME.out.transcript_qc_stats
                .map{ _meta, stats -> [stats]}
                .collect()
                .ifEmpty([])
        )
    }

    // oarfish expects deduplicated reads
    if (transcript_quants) {
        PROCESS_LONGREAD_SCRNA_TRANSCRIPT (
            transcript_fasta,
            transcript_fai,
            gtf,
            ch_extracted_fastq_cdna,
            ch_rseqc_bed,
            ch_corrected_bc_info_cdna,
            transcript_quants,
            params.dedup_tool,
            false, // Indicates this is NOT genome aligned
            params.fasta_delimiter,
            params.skip_save_minimap2_index,
            params.skip_qc,
            true, // RSeQC does not work well with transcriptome alignments
            true, // Nanocomp does not work well with transcriptome alignments
            params.skip_seurat,
            false // Oarfish requires deduplication, so cannot skip it
        )

        ch_versions = ch_versions.mix(PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.versions)

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_flagstat.collect{it -> it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_rseqc_read_dist.collect{it -> it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.minimap_nanocomp_bam_txt.collect{it -> it[1]}.ifEmpty([])
        )

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.bc_tagged_flagstat.collect{it -> it[1]}.ifEmpty([])
        )

        // The transcriptome arm always deduplicates (see above), so its dedup QC is
        // always reported.
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.dedup_flagstat.collect{it -> it[1]}.ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.dedup_summary.collect{it -> it[1]}.ifEmpty([])
        )

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            PROCESS_LONGREAD_SCRNA_TRANSCRIPT.out.transcript_qc_stats
                .map{ _meta, stats -> [stats]}
                .collect()
                .ifEmpty([])
        )
    }

    //
    // SUBWORKFLOW: Quantify cDNA with lr-kallisto, straight from the fastq
    //

    if (fastq_quants) {
        QUANTIFY_SCRNA_LRKALLISTO (
            ch_extracted_fastq_cdna,
            ch_corrected_bc_info_cdna,
            ch_barcode_counts_cdna,
            PREPARE_REFERENCE_FILES.out.derived_transcript_fasta,
            genome_fasta,
            gtf,
            params.skip_qc,
            params.skip_seurat
        )

        ch_versions = ch_versions.mix(QUANTIFY_SCRNA_LRKALLISTO.out.versions)

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            QUANTIFY_SCRNA_LRKALLISTO.out.gene_qc_stats
                .map{ _meta, stats -> [stats]}
                .collect()
                .ifEmpty([])
        )
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
            QUANTIFY_SCRNA_LRKALLISTO.out.transcript_qc_stats
                .map{ _meta, stats -> [stats]}
                .collect()
                .ifEmpty([])
        )
    }

    //
    // SUBWORKFLOW: Align and deduplicate DNA samples
    //

    ALIGN_DEDUPLICATE_DNA (
        genome_fasta,
        genome_fai,
        ch_extracted_fastq_dna,
        params.skip_save_minimap2_index,
        params.skip_qc,
        params.skip_bam_nanocomp,
        params.skip_dedup
    )

    ch_versions = ch_versions.mix(ALIGN_DEDUPLICATE_DNA.out.versions)

    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ALIGN_DEDUPLICATE_DNA.out.flagstat.collect{it[1]}.ifEmpty([])
    )
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ALIGN_DEDUPLICATE_DNA.out.stats.collect{it -> it[1]}.ifEmpty([])
    )
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ALIGN_DEDUPLICATE_DNA.out.idxstats.collect{it -> it[1]}.ifEmpty([])
    )
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ALIGN_DEDUPLICATE_DNA.out.dedup_metrics.collect{it[1]}.ifEmpty([])
    )

    //
    // The per-sample *.dedup_summary.tsv is long (metric<TAB>value), which MultiQC
    // custom content cannot transpose. Collect one wide row per sample instead: a
    // table of the main counts and rates, and the duplicate rates for the general
    // statistics. MultiQC reads each file into one section only, hence two files.
    // The rows are named <id>.dna, as the samtools stats of the same bam are.
    //
    def ch_dna_dedup_rows = ALIGN_DEDUPLICATE_DNA.out.dedup_summary
        .map { meta, summary -> [ "${meta.id}.${meta.type}", dnaDedupMetrics(summary) ] }
        .filter { _sample, metrics -> metrics }

    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ch_dna_dedup_rows
            .map { sample, metrics -> wideRow(sample, metrics, dnaDedupTableColumns()) }
            .collectFile(name: 'dna_dedup_summary.tsv', seed: wideHeader(dnaDedupTableColumns()), sort: true)
    )
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ch_dna_dedup_rows
            .map { sample, metrics -> wideRow(sample, metrics, dnaDedupGeneralStatsColumns()) }
            .collectFile(name: 'dna_dedup_general_stats.tsv', seed: wideHeader(dnaDedupGeneralStatsColumns()), sort: true)
    )
    ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(
        ALIGN_DEDUPLICATE_DNA.out.nanocomp_bam_txt.collect{it[1]}.ifEmpty([])
    )

    //
    // SOFTWARE_VERSIONS
    //

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'lrsinglecell_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )

    ch_multiqc_list = channel.empty()
    if (!params.skip_qc && !params.skip_multiqc){

        //
        // MultiQC configs, in the order MultiQC applies them (a later one wins): the
        // pipeline's, the report comment for this version, then --multiqc_config.
        //
        def ch_multiqc_report_comment = channel.of(multiqcReportComment())
            .collectFile(name: 'multiqc_report_comment.yml')
        def ch_multiqc_configs = ch_multiqc_config
            .concat(ch_multiqc_report_comment, ch_multiqc_custom_config)
            .collect()
            .map { cfg -> [ cfg ] }

        // Raw-read QC runs whenever any of its inputs does: NanoPlot (on by
        // default), FastQC or the fastq NanoComp.
        if (!params.skip_fastqc || !params.skip_fastq_nanocomp || !params.skip_nanoplot) {
            //
            // MODULE: MultiQC for raw data
            //
            ch_multiqc_rawqc_files = ch_fastqc_multiqc_pretrim.collect().ifEmpty([])
            // ch_nanostat_pretrim is ifEmpty([]): with NanoPlot skipped it emits a
            // bare [], and it[1] on that would hand MultiQC a null path.
            if (!params.skip_nanoplot) {
                ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanostat_pretrim.collect{it -> it[1]}.ifEmpty([]))
            }
            ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanocomp_fastq_txt.collect{it -> it[1]}.ifEmpty([]))
            ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_collated_versions.collect())

            MULTIQC_RAWQC (
                channel.value([ [ 'id': 'rawqc' ] ])
                    .combine(ch_multiqc_rawqc_files.collect().map { files -> [ files ] })
                    .combine(ch_multiqc_configs)
                    .combine(ch_multiqc_logo.collect().ifEmpty([[]]).map { logo -> [ logo ] })
                    .combine(channel.value([ [] ]))
                    .combine(channel.value([ [] ]))
                    .map { meta, files, config, logo, replace, sample ->
                        [ meta, files, config, logo.flatten(), replace, sample ]
                    }
            )
        }

        //
        // MODULE: MultiQC for final pipeline outputs
        //
        def summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

        // The methods text names the tools of the arms this samplesheet runs
        def ch_methods_description = ch_samplesheet
            .map { meta, _fastqs -> meta.type }
            .unique()
            .toSortedList()
            .map { sample_types -> [ sample_types ] }
            .combine(ch_multiqc_custom_methods_description)
            .map { sample_types, methods_yaml -> methodsDescriptionText(methods_yaml, sample_types) }

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_collated_versions)
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postextract.collect().ifEmpty([]))

        MULTIQC_FINALQC (
            channel.value([ [ 'id': 'finalqc' ] ])
                .combine(ch_multiqc_finalqc_files.collect().map { files -> [ files ] })
                .combine(ch_multiqc_configs)
                .combine(ch_multiqc_logo.collect().ifEmpty([[]]).map { logo -> [ logo ] })
                .combine(channel.value([ [] ]))
                .combine(channel.value([ [] ]))
                .map { meta, files, config, logo, replace, sample ->
                    [ meta, files, config, logo.flatten(), replace, sample ]
                }
        )
        ch_multiqc_report = MULTIQC_FINALQC.out.report
        ch_multiqc_list = ch_multiqc_report.map { _meta, report -> [report] }.toList()
    }

    emit:
    multiqc_report = ch_multiqc_list
    versions       = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Rows of mark_dna_duplicates.py's *.dedup_summary.tsv shown in the MultiQC DNA
// deduplication table, in the order the script writes them. The headers in
// assets/multiqc_config.yml describe each one. A row missing from the summary is
// left empty, so the table does not depend on every row being there.
//
def dnaDedupTableColumns() {
    return [
        'primary_reads',
        'molecules',
        'duplicates',
        'duplicate_rate',
        'nuclear_duplicate_rate',
        'exact_only_contigs_primary_reads',
        'exact_only_contigs_duplicates',
        'dup_by_key_5prime_strand',
        'dup_added_by_same_raw_read',
        'dup_added_by_5prime_jitter',
        'dup_added_by_aa_opposite_end',
        'no_barcode',
    ]
}

// The rates shown in the MultiQC general statistics
def dnaDedupGeneralStatsColumns() {
    return [ 'nuclear_duplicate_rate', 'duplicate_rate' ]
}

// The metric<TAB>value rows of a *.dedup_summary.tsv as a map; empty for a stub file
def dnaDedupMetrics(summary) {
    def metrics = [:]
    summary.readLines().drop(1).each { line ->
        def fields = line.split('\t', -1)
        if (fields.size() == 2) {
            metrics[fields[0]] = fields[1]
        }
    }
    return metrics
}

def wideHeader(columns) {
    return ([ 'sample' ] + columns).join('\t') + '\n'
}

def wideRow(sample, metrics, columns) {
    return ([ sample ] + columns.collect { column -> metrics.get(column, '') }).join('\t') + '\n'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Performs feature quantification for long read single-cell rna data
//

include { SAMTOOLS_FILTER_SORT } from '../../../modules/local/samtools_filter_sort'
include { OARFISH              } from '../../../modules/local/oarfish'
include { QC_SCRNA             } from '../../../subworkflows/local/qc_scrna'

workflow QUANTIFY_SCRNA_OARFISH {
    take:
        ch_bam      // channel: [ val(meta), path(bam) ]
        ch_flagstat // channel: [ val(meta), path(flagstat) ]

        val_skip_qc     // bool: Skip qc steps
        val_skip_seurat // bool: Skip seurat qc steps

    main:
        ch_versions = channel.empty()

        //
        // MODULE: Drop reads without a called-cell barcode and sort by barcode
        //
        // flexiplex -a true leaves CB as '-' on reads that matched no known
        // barcode. oarfish groups on CB, so without this filter all of them
        // would pool into a single '-' pseudo-cell. Filtering in the sort that
        // oarfish needs anyway costs no extra pass over the BAM.
        //
        SAMTOOLS_FILTER_SORT ( ch_bam )
        ch_versions = ch_versions.mix(SAMTOOLS_FILTER_SORT.out.versions_samtools_filter_sort)

        //
        // MODULE: Oarfish
        //
        OARFISH ( SAMTOOLS_FILTER_SORT.out.bam )
        ch_versions = ch_versions.mix(OARFISH.out.versions_oarfish)

        ch_transcript_qc_stats = channel.empty()
        if (!val_skip_qc && !val_skip_seurat) {
            QC_SCRNA(
                OARFISH.out.features
                    .join(OARFISH.out.barcodes, by: [0])
                    .join(OARFISH.out.mtx, by: [0])
                    .map{
                        meta,features,barcodes,mtx ->
                            def new_meta = meta
                            [ new_meta, [ features, barcodes, mtx ]]
                    },
                ch_flagstat,
                "MEX"
            )
            ch_transcript_qc_stats = QC_SCRNA.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA.out.versions)
        }

    emit:
        versions            = ch_versions            // channel: [ val(meta), path(versions) ]
        transcript_qc_stats = ch_transcript_qc_stats // channel: [ val(meta), path(transcript_qc_stats) ]
}

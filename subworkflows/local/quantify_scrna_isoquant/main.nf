//
// Performs feature quantification for long read single-cell rna data
//

include { ISOQUANT                                } from '../../../modules/local/isoquant'
include { DROP_BARCODE as DROP_BARCODE_GENE       } from '../../../modules/local/drop_barcode'
include { DROP_BARCODE as DROP_BARCODE_TRANSCRIPT } from '../../../modules/local/drop_barcode'
include { QC_SCRNA as QC_SCRNA_GENE               } from '../../../subworkflows/local/qc_scrna'
include { QC_SCRNA as QC_SCRNA_TRANSCRIPT         } from '../../../subworkflows/local/qc_scrna'

workflow QUANTIFY_SCRNA_ISOQUANT {
    take:
        ch_bam      // channel: [ val(meta), path(bam) ]
        ch_bai      // channel: [ val(meta), path(bai) ]
        ch_flagstat // channel: [ val(meta), path(flagstat) ]
        ch_fasta    // channel: [ val(meta), path(fasta) ]
        ch_fai      // channel: [ val(meta), path(fai) ]
        ch_gtf      // channel: [ val(meta), path(gtf) ]

        val_group_category // str: Isoquant --read_group, e.g. 'tag:CB' or 'tag:XB'
        val_drop_barcode   // str: barcode to remove from the matrices (e.g. '-'), or '' to keep them as IsoQuant wrote them
        val_skip_qc        // bool: Skip qc steps
        val_skip_seurat    // bool: Skip seurat qc steps

    main:
        ch_versions = channel.empty()

        isoquant_input = ch_bam
            .join(ch_bai, by: [0])
            .combine(ch_fasta.map {_meta, fasta -> [fasta]})
            .combine(ch_fai.map {_meta, fai -> [fai]})
            .combine(ch_gtf.map {_meta, gtf -> [gtf]})

        //
        // MODULE: Isoquant
        //
        ISOQUANT (
            isoquant_input,
            val_group_category
        )
        ch_versions = ch_versions.mix(ISOQUANT.out.versions_isoquant)

        //
        // MODULE: Drop a barcode, e.g. flexiplex's '-' placeholder for reads that matched
        // no known barcode, which a CB-grouped run otherwise reports as one pseudo-cell.
        // It runs before publishing and Seurat QC, so neither sees that column.
        //
        ch_gene_mtx                = ISOQUANT.out.grouped_gene_mtx
        ch_gene_mtx_barcodes       = ISOQUANT.out.grouped_gene_mtx_barcodes
        ch_gene_mtx_features       = ISOQUANT.out.grouped_gene_mtx_features
        ch_transcript_mtx          = ISOQUANT.out.grouped_transcript_mtx
        ch_transcript_mtx_barcodes = ISOQUANT.out.grouped_transcript_mtx_barcodes
        ch_transcript_mtx_features = ISOQUANT.out.grouped_transcript_mtx_features
        if (val_drop_barcode) {
            DROP_BARCODE_GENE (
                ISOQUANT.out.grouped_gene_mtx
                    .join(ISOQUANT.out.grouped_gene_mtx_barcodes, by: [0])
                    .join(ISOQUANT.out.grouped_gene_mtx_features, by: [0]),
                val_drop_barcode
            )
            ch_gene_mtx          = DROP_BARCODE_GENE.out.mtx
            ch_gene_mtx_barcodes = DROP_BARCODE_GENE.out.barcodes
            ch_gene_mtx_features = DROP_BARCODE_GENE.out.features
            ch_versions = ch_versions.mix(DROP_BARCODE_GENE.out.versions_drop_barcode)

            DROP_BARCODE_TRANSCRIPT (
                ISOQUANT.out.grouped_transcript_mtx
                    .join(ISOQUANT.out.grouped_transcript_mtx_barcodes, by: [0])
                    .join(ISOQUANT.out.grouped_transcript_mtx_features, by: [0]),
                val_drop_barcode
            )
            ch_transcript_mtx          = DROP_BARCODE_TRANSCRIPT.out.mtx
            ch_transcript_mtx_barcodes = DROP_BARCODE_TRANSCRIPT.out.barcodes
            ch_transcript_mtx_features = DROP_BARCODE_TRANSCRIPT.out.features
            ch_versions = ch_versions.mix(DROP_BARCODE_TRANSCRIPT.out.versions_drop_barcode)
        }

        ch_gene_qc_stats = channel.empty()
        ch_transcript_qc_stats = channel.empty()
        if (!val_skip_qc && !val_skip_seurat){
            QC_SCRNA_GENE (
                ch_gene_mtx_features
                    .join(ch_gene_mtx_barcodes, by: [0])
                    .join(ch_gene_mtx, by: [0])
                    // Keep meta unchanged: QC_SCRNA joins on it with ch_flagstat,
                    // whose meta still carries type/cell_count/single_end.
                    .map{
                        meta,features,barcodes,mtx ->
                            [ meta, [ features, barcodes, mtx ]]
                    },
                ch_flagstat,
                "MEX"
            )
            ch_gene_qc_stats = QC_SCRNA_GENE.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_GENE.out.versions)

            QC_SCRNA_TRANSCRIPT (
                ch_transcript_mtx_features
                    .join(ch_transcript_mtx_barcodes, by: [0])
                    .join(ch_transcript_mtx, by: [0])
                    // Keep meta unchanged: QC_SCRNA joins on it with ch_flagstat,
                    // whose meta still carries type/cell_count/single_end.
                    .map{
                        meta,features,barcodes,mtx ->
                            [ meta, [ features, barcodes, mtx ]]
                    },
                ch_flagstat,
                "MEX"
            )
            ch_transcript_qc_stats = QC_SCRNA_TRANSCRIPT.out.seurat_stats
            ch_versions = ch_versions.mix(QC_SCRNA_TRANSCRIPT.out.versions)
        }

    emit:
        versions                 = ch_versions                // channel: [ val(meta), path(versions) ]
        gene_features_file       = ch_gene_mtx_features       // channel: [ val(meta), path(gene_features) ]
        gene_barcodes_file       = ch_gene_mtx_barcodes       // channel: [ val(meta), path(gene_barcodes) ]
        gene_mtx_file            = ch_gene_mtx                // channel: [ val(meta), path(gene_mtx) ]
        transcript_features_file = ch_transcript_mtx_features // channel: [ val(meta), path(transcript_features) ]
        transcript_barcodes_file = ch_transcript_mtx_barcodes // channel: [ val(meta), path(transcript_barcodes) ]
        transcript_mtx_file      = ch_transcript_mtx          // channel: [ val(meta), path(transcript_mtx) ]
        gene_qc_stats            = ch_gene_qc_stats           // channel: [ val(meta), path(gene_qc_stats) ]
        transcript_qc_stats      = ch_transcript_qc_stats     // channel: [ val(meta), path(transcript_qc_stats) ]
}

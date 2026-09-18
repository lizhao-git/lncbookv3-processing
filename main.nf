nextflow.enable.dsl = 2

include { REFERENCE } from './subworkflows/local/reference/main'
include { CLINVAR } from './subworkflows/local/clinvar/main'
include { COSMIC } from './subworkflows/local/cosmic/main'
include { GWAS } from './subworkflows/local/gwas/main'
include { SMPROT } from './subworkflows/local/smprot/main'
include { METHYLATION } from './subworkflows/local/methylation/main'
include { CERNA } from './subworkflows/local/cerna/main'
include { CONSERVATION } from './subworkflows/local/conservation/main'

workflow LNCBOOKV3_PROCESSING {
    take:
    gtf_ch

    main:
    REFERENCE(gtf_ch)

    clinvar_annotations_ch = Channel.empty()
    cosmic_annotations_ch = Channel.empty()
    gwas_annotations_ch = Channel.empty()
    smprot_annotations_ch = Channel.empty()
    methylation_outputs_ch = Channel.empty()
    cerna_outputs_ch = Channel.empty()
    conservation_outputs_ch = Channel.empty()

    if (params.run_clinvar && params.clinvar_vcf) {
        clinvar_ch = Channel.fromPath(params.clinvar_vcf, checkIfExists: true).map { file -> tuple([id: file.baseName], file) }
        // --annotation accepts gtf/gff3/bed (auto-detected); falls back to --gtf
        def annotation_path = params.annotation ?: params.gtf
        annotation_ch = Channel.fromPath(annotation_path, checkIfExists: true).map { file -> tuple([id: file.baseName], file) }
        CLINVAR(clinvar_ch, annotation_ch)
        clinvar_annotations_ch = CLINVAR.out.annotations
    }

    if (params.run_cosmic && params.cosmic_tsv) {
        cosmic_ch = Channel.fromPath(params.cosmic_tsv, checkIfExists: true).map { file -> tuple([id: file.baseName], file) }
        // --annotation accepts gtf/gff3/bed (auto-detected); falls back to --gtf
        def cosmic_annotation_path = params.annotation ?: params.gtf
        cosmic_annotation_ch = Channel.fromPath(cosmic_annotation_path, checkIfExists: true).map { file -> tuple([id: file.baseName], file) }
        COSMIC(cosmic_ch, cosmic_annotation_ch)
        cosmic_annotations_ch = COSMIC.out.annotations
    }

    if (params.run_gwas_catalog && params.gwas_tsv) {
        gwas_ch = Channel.fromPath(params.gwas_tsv, checkIfExists: true).map { file -> tuple([id: file.baseName], file) }
        // --annotation accepts gtf/gff3/bed (auto-detected); falls back to --gtf
        def gwas_annotation_path = params.annotation ?: params.gtf
        gwas_annotation_ch = Channel.fromPath(gwas_annotation_path, checkIfExists: true).map { file -> tuple([id: file.baseName], file) }
        GWAS(gwas_ch, gwas_annotation_ch)
        gwas_annotations_ch = GWAS.out.annotations
    }

    if (params.run_smprot && params.smprot_tsv) {
        smprot_files = params.smprot_tsv instanceof List ? params.smprot_tsv.collect { file(it) } : [file(params.smprot_tsv)]
        smprot_ch = Channel.of(tuple([id: 'smprot'], smprot_files))
        // --annotation accepts gtf/gff3/bed (auto-detected); falls back to --gtf
        def smprot_annotation_path = params.annotation ?: params.gtf
        smprot_annotation_ch = Channel.fromPath(smprot_annotation_path, checkIfExists: true).map { file -> tuple([id: file.baseName], file) }
        SMPROT(smprot_ch, smprot_annotation_ch)
        smprot_annotations_ch = SMPROT.out.annotations
    }

    if (params.run_methylation && params.methylation_manifest) {
        methylation_manifest_ch = Channel.of(tuple([id: 'methylation'], file(params.methylation_manifest)))
        methylation_root_ch = Channel.value(params.methylation_data_root ? file(params.methylation_data_root) : [])
        METHYLATION(methylation_manifest_ch, methylation_root_ch)
        methylation_outputs_ch = METHYLATION.out.outputs_dir
    }

    if (params.run_cerna && params.cerna_manifest) {
        cerna_manifest_ch = Channel.of(tuple([id: 'cerna'], file(params.cerna_manifest)))
        cerna_root_ch = Channel.value(params.cerna_data_root ? file(params.cerna_data_root) : [])
        CERNA(cerna_manifest_ch, cerna_root_ch)
        cerna_outputs_ch = CERNA.out.outputs_dir
    }

    if (params.run_conservation && params.conservation_manifest) {
        conservation_manifest_ch = Channel.of(tuple([id: 'conservation'], file(params.conservation_manifest)))
        conservation_root_ch = Channel.value(params.conservation_data_root ? file(params.conservation_data_root) : [])
        CONSERVATION(conservation_manifest_ch, conservation_root_ch)
        conservation_outputs_ch = CONSERVATION.out.outputs_dir
    }

    emit:
    validated_gtf = REFERENCE.out.validated_gtf
    gtf_validation_report = REFERENCE.out.gtf_validation_report
    gtf_features = REFERENCE.out.features
    gtf_features_validation_report = REFERENCE.out.features_validation_report
    clinvar_annotations = clinvar_annotations_ch
    cosmic_annotations = cosmic_annotations_ch
    gwas_annotations = gwas_annotations_ch
    smprot_annotations = smprot_annotations_ch
    methylation_outputs = methylation_outputs_ch
    cerna_outputs = cerna_outputs_ch
    conservation_outputs = conservation_outputs_ch
}

workflow {
    if (!params.gtf) {
        error "Missing required parameter: --gtf"
    }
    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true).map { file -> tuple([id: 'reference'], file) }
    LNCBOOKV3_PROCESSING(gtf_ch)
}

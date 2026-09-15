nextflow.enable.dsl = 2

include { METHYLATION as RUN_METHYLATION } from '../../../modules/local/methylation/main'

workflow METHYLATION {
    take:
    manifest_ch
    data_root_ch

    main:
    RUN_METHYLATION(manifest_ch, data_root_ch)

    emit:
    outputs_dir = RUN_METHYLATION.out.outputs_dir
    report = RUN_METHYLATION.out.report
}

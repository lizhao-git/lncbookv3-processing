nextflow.enable.dsl = 2

include { CONSERVATION as RUN_CONSERVATION } from '../../../modules/local/conservation/main'

workflow CONSERVATION {
    take:
    manifest_ch
    data_root_ch

    main:
    RUN_CONSERVATION(manifest_ch, data_root_ch)

    emit:
    outputs_dir = RUN_CONSERVATION.out.outputs_dir
    report = RUN_CONSERVATION.out.report
}

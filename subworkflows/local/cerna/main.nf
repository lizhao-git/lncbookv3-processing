nextflow.enable.dsl = 2

include { CERNA as RUN_CERNA } from '../../../modules/local/cerna/main'

workflow CERNA {
    take:
    manifest_ch
    data_root_ch

    main:
    RUN_CERNA(manifest_ch, data_root_ch)

    emit:
    outputs_dir = RUN_CERNA.out.outputs_dir
}

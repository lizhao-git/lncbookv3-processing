nextflow.enable.dsl = 2

include { CERNA_PREDICT } from '../../../modules/local/cerna_predict/main'

workflow CERNA_PREDICT {
    take:
    manifest_ch
    data_root_ch

    main:
    CERNA_PREDICT(manifest_ch, data_root_ch)

    emit:
    outputs_dir = CERNA_PREDICT.out.outputs_dir
}

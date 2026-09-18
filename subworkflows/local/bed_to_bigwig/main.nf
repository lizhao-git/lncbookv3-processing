nextflow.enable.dsl = 2

// BED -> bigWig subworkflow
//
// Chains bed_to_bedgraph (bedtools genomecov -bg + sort) then
// ucsc/bedgraphtobigwig to produce a bigWig coverage track from BED
// intervals in one self-contained call.
//
// Usage:
//   include { BED_TO_BIGWIG } from './subworkflows/local/bed_to_bigwig/main'
//
//   workflow {
//       bed_ch = Channel.of(tuple([id: 'sample'], file('input.bed')))
//       sizes  = file('chrom.sizes')
//       BED_TO_BIGWIG(bed_ch, sizes)
//   }

include { BED_TO_BEDGRAPH as RUN_BED_TO_BEDGRAPH } from '../../../modules/local/bed_to_bedgraph/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BED_TO_BIGWIG {
    take:
    bed_ch       // channel: tuple val(meta), path(bed)
    chrom_sizes  // channel: path(chrom.sizes)

    main:
    // Generate sorted bedGraph coverage from BED intervals
    RUN_BED_TO_BEDGRAPH(bed_ch, chrom_sizes, 'coverage.bedgraph')

    // Convert sorted bedGraph to bigWig
    UCSC_BEDGRAPHTOBIGWIG(RUN_BED_TO_BEDGRAPH.out.bedgraph, chrom_sizes)

    emit:
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    bedgraph = RUN_BED_TO_BEDGRAPH.out.bedgraph
}

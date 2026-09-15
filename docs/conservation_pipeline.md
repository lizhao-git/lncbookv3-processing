# lncRNA Conservation Analysis Pipeline

This document organizes the lncRNA conservation workflow from the supplied
figure and the two post-processing scripts:

- `/Volumes/Extreme SSD/Projects/group_process_filter_lnc.R`
- `/Volumes/Extreme SSD/Projects/group_process_filter_pc.R`

The pipeline follows the same conceptual design used by LncBook 2.0: identify
cross-species transcript mappings from UCSC genome assemblies, annotations and
chain files; integrate sequence similarity, transcript coverage, exon/intron
structure, genomic context and coding potential; then derive sequence
conservation, homologous relationships, synteny and gene age.

## Reference Method

LncBook 2.0 used pairwise genome mappings from human to 40 vertebrate species
based on UCSC resources. Its conservation analysis is built around three
quantitative filters:

- mapped alignment length greater than or equal to 50 bp;
- transcript coverage greater than or equal to 0.2;
- intron coverage greater than or equal to the Q50 threshold.

The homologous gene table in LncBook 2.0 reports each human lncRNA gene and
its orthologous or homologous genes in other species. Gene age is assigned from
the most ancient species in which a homolog is detected.

## High-Level Flow

```text
genome sequence + chain file + annotation
  |
  |-- prepare human/query transcript BED and metadata
  |-- prepare target species annotation metadata
  |
  |-- pslMap / liftover-style transcript mapping
        |
        |-- query.bed + target.bed
        |-- query.fa + target.fa
        |
        |-- segment-level sequence identity
        |-- duplicate hit filtering and grouping
        |-- target transcript/exon/CDS/ORF annotation merge
        |
        |-- two-way information integration
              |
              |-- sequence conservation
              |-- gene homology
              |-- synteny
              |-- gene age
```

## Stage Design

### A. Input Preparation

Inputs:

- Human lncRNA annotation: GTF/GFF3/BED12.
- Human genome FASTA.
- Target species genome FASTA.
- Human-to-target UCSC chain file.
- Target species annotation: GTF/GFF3/BED12.
- Optional ORF/CDS annotation for protein-coding controls and target
  transcript context.

Outputs:

- `query_transcripts.bed`: human lncRNA transcript intervals in BED12.
- `query_exons.bed`: human lncRNA exon intervals.
- `query_lengths.tsv`: transcript/exon block lengths used for coverage.
- `query_meta.tsv`: transcript/gene metadata.
- `target_transcripts.bed`: target species transcript intervals.
- `target_exons.bed`: target species exon intervals.
- `target_cds_or_orf.tsv`: target coding or ORF context.
- `target_meta.tsv`: target transcript metadata.

Recommended implementation:

- Reuse `scripts/genomic_intervals/annotation.py` for GTF/GFF3 validation and
  BED-like interval extraction.
- Add a conservation-specific converter that emits BED12, exon positions,
  transcript length and metadata files in the same schema expected by the R
  scripts.

### B. Cross-Species Mapping

Inputs:

- `query_transcripts.bed`
- `target_transcripts.bed`
- Human and target genome FASTA
- Chain file

Operation:

- Run `pslMap` or the equivalent UCSC chain-based mapping to project human
  transcript intervals to target genome coordinates.
- Extract paired query and target mapped segment FASTA.

Outputs:

- `{species}_lnc_query_simply.fa`
- `{species}_lnc_target_simply.fa`
- `{species}_lnc_group_info.txt`
- `{species}_lnc_lengths.txt`
- `{species}_lnc_all_target_info.txt`
- `{species}_exoncount_ORFexonrange.txt`

Notes:

- The figure emphasizes this as the central mapping stage.
- The paired FASTA identifiers must remain one-to-one after sorting; both R
  scripts explicitly check name and length equality.

### C. Segment Similarity and Hit Filtering

Implemented by:

- `group_process_filter_lnc.R` for lncRNA mappings.
- `group_process_filter_pc.R` for protein-coding control mappings.

Inputs consumed by the lncRNA script:

- `resourse/fas/{prefix}_lnc_query_simply.fa`
- `resourse/fas/{prefix}_lnc_target_simply.fa`
- `resourse/lnc_group_info/{prefix}_lnc_group_info.txt`
- `resourse/lnc_lengths/{prefix}_lnc_lengths.txt`
- `resourse/bedtool_exonCDS/{prefix}_lnc_all_target_info.txt`
- `resourse/CDSexon_range/{prefix}_exoncount_ORFexonrange.txt`
- `resourse/exon_positions.txt`
- `resourse/all_lnc_length.txt`
- `resourse/hg38_all_trans_meta_info_new.txt`

Core operations:

- Sort query and target mapped FASTA records by sequence ID.
- Check one-to-one query-target ID mapping.
- Check mapped sequence lengths are identical.
- Compute per-segment `matched_length` and `match_ratio`.
- Collapse duplicated segment records back to transcript/group IDs.
- Remove overlapping duplicate hits when two hits overlap by at least 50% of
  the smaller mapped span, retaining the hit with greater matched length.
- Regroup nearby duplicated hits when they are on the same strand and within
  100 kb.
- Merge target transcript/exon/ORF metadata.
- Compute exon coverage:
  - `exon_counts`
  - `exon_covered`
  - `exon_covered_num`
  - `exon_covered_ratio`
- Compute true non-redundant mapped length across grouped segments.
- Merge transcript length and human transcript metadata.

Output:

- `{prefix}_lnc_all_statistics_new.txt`

Important output fields:

- `ids`: original human transcript ID.
- `grouped_ids`: grouped conserved segment ID.
- `transcript_length`: original transcript length.
- `true_mapped_length`: non-redundant mapped length.
- `mapped_length`: total mapped segment length.
- `matched_length`: identical bases across mapped segments.
- `match_ratio`: `matched_length / mapped_length`.
- `exon_counts`, `exon_covered`, `exon_covered_num`,
  `exon_covered_ratio`.
- target transcript/exon/ORF context columns.
- human transcript/gene metadata columns.

### D. Protein-Coding Control Processing

Implemented by:

- `group_process_filter_pc.R`

Inputs are parallel to the lncRNA script, with `pc` resource directories:

- `resourse/fas/{prefix}_pc_query_simply.fa`
- `resourse/fas/{prefix}_pc_target_simply.fa`
- `resourse/pc_group_info/{prefix}_pc_group_info.txt`
- `resourse/pc_lengths/{prefix}_pc_lengths.txt`
- `resourse/modified_pc_CDS_positions/{prefix}_pc_query_add_CDS_positions_simply_extractedCDSpositions.txt`
- `resourse/Final_pc_query_beds/{prefix}_pc_CDS_start_end_cover_index.txt`
- `resourse/bedtool_exonCDS/{prefix}_all_target_info.txt`
- `resourse/CDSexon_range/{prefix}_exoncount_ORFexonrange.txt`
- `resourse/exon_positions_pro.txt`
- `resourse/all_pro_length.txt`
- `resourse/all_pro_CDS_length.txt`
- `resourse/Final_pc_query_beds/LncBook_pc_CDS_positions.txt`
- `resourse/hg38_all_trans_meta_info_new.txt`

Additional protein-coding operations:

- Compute mapped and matched length inside CDS.
- Compute `CDS_match_ratio`.
- Merge CDS start/end cover index.
- Compute total CDS length with phase correction.
- Compute `CDS_coverage`.

Output:

- `{prefix}_pc_all_statistics_new.txt`

The PC branch is useful as a positive control and for calibrating conservation
thresholds, especially CDS-aware metrics that are not applicable to lncRNAs.

### E. Conservation Classification

Inputs:

- `{species}_lnc_all_statistics_new.txt`
- Species tree/order table.
- Per-species Q50 intron coverage threshold.

Recommended filters:

- `true_mapped_length >= 50`
- `transcript_coverage = true_mapped_length / transcript_length >= 0.2`
- intron coverage greater than or equal to the species-specific or global Q50
  threshold.

Outputs:

- `lncRNA_sequence_conservation.tsv`
- `lncRNA_conservation_segments.bed`
- `lncRNA_conservation_summary.tsv`

Recommended fields:

- `human_gene_id`
- `human_transcript_id`
- `target_species`
- `target_gene_id`
- `target_transcript_id`
- `mapped_length`
- `true_mapped_length`
- `matched_length`
- `match_ratio`
- `transcript_length`
- `transcript_coverage`
- `exon_covered_ratio`
- `intron_coverage`
- `pass_length_filter`
- `pass_transcript_coverage_filter`
- `pass_intron_coverage_filter`
- `is_conserved`

### F. Gene Homology

Inputs:

- Conserved lncRNA transcript/species pairs from stage E.
- Human transcript-to-gene mapping.
- Target transcript-to-gene mapping.

Operation:

- Collapse transcript-level mappings to gene-level homologous relationships.
- Prefer mappings with stronger evidence when multiple target genes exist:
  higher `true_mapped_length`, higher `match_ratio`, higher transcript
  coverage, and broader exon coverage.

Outputs:

- `lncRNA_gene_homology.tsv`
- `lncRNA_gene_homology_best_hit.tsv`

### G. Synteny

Inputs:

- Gene homology table.
- Human and target neighboring gene annotations.
- Chain-derived genomic mapping or mapped transcript loci.

Operation:

- For each conserved lncRNA locus, compare neighboring protein-coding genes or
  stable flanking anchors between human and target species.
- Mark syntenic support when mapped lncRNA homologs preserve local gene order
  or neighborhood anchors.

Outputs:

- `lncRNA_synteny.tsv`
- `lncRNA_syntenic_homology.tsv`

### H. Gene Age

Inputs:

- Gene homology table.
- Ordered species phylogeny or clade table.

Operation:

- For each human lncRNA gene, find the most ancient species or clade where at
  least one homolog is detected.
- Assign this species/clade as the gene age.

Outputs:

- `lncRNA_gene_age.tsv`
- `lncRNA_gene_age_summary.tsv`

Recommended fields:

- `human_gene_id`
- `human_gene_name`
- `oldest_detected_species`
- `oldest_detected_clade`
- `supporting_target_genes`
- `supporting_transcripts`
- `species_count`

## How the Supplied R Scripts Fit

The two supplied R scripts are not the full conservation pipeline. They are
post-pslMap processing steps.

They assume upstream files have already been created:

- mapped query and target FASTA;
- mapped BED/group interval files;
- original transcript length files;
- exon position files;
- target transcript/CDS/ORF annotation files;
- human transcript metadata.

They produce the per-species, per-prefix statistics table used by downstream
conservation classification. The lncRNA script is the primary branch for
LncBook lncRNA conservation; the protein-coding script is a comparable control
branch with additional CDS metrics.

## Recommended Modularization for This Repository

The repository now contains a new `scripts/conservation/` package with these
stages:

- `postprocess.py`: Python replacement for the sequence identity, duplicate
  removal, same-strand 100 kb grouping, exon coverage and context merge logic
  in `group_process_filter_lnc.R` and `group_process_filter_pc.R`. BED12 block
  parsing and interval arithmetic are reused from the shared
  `scripts/format_convert/genomic_intervals/bed12.py` layer.
- `group_process_filter_lnc.py`: compatibility wrapper for the lncRNA branch.
- `group_process_filter_pc.py`: compatibility wrapper for the protein-coding
  control branch.
- `classify.py`: applies LncBook 2.0 style length, transcript coverage and
  intron coverage filters.
- `homology.py`: collapses transcript-level evidence to gene-level homology
  and best hits.
- `synteny.py`: annotates homologs with optional synteny anchor support.
- `gene_age.py`: assigns oldest detected species/clade from a species order
  table.
- `pipeline.py`: manifest-driven orchestrator used by Nextflow.

Current Nextflow modules:

- `modules/local/conservation/main.nf`: runs the manifest-driven conservation pipeline.
- `subworkflows/local/conservation/main.nf`: standalone conservation branch workflow.
- `modules/local/pslmap_chain/main.nf`: optional upstream wrapper for UCSC `pslMap`
  when PSL mappings and chain files are generated in a Nextflow run.

The species-level mapping/post-processing workflow should scatter over a
species manifest. Each manifest row should contain:

- `species`
- `assembly`
- `genome_fasta`
- `annotation`
- `annotation_format`
- `chain_from_human`
- optional `target_orf_or_cds`
- optional species-specific intron Q50 threshold.

## Minimal End-to-End Outputs

For database import, the pipeline should at minimum produce:

- `conservation/segments/{species}_lnc_segments.tsv`
- `conservation/statistics/{species}_lnc_all_statistics.tsv`
- `conservation/sequence_conservation.tsv`
- `conservation/gene_homology.tsv`
- `conservation/synteny.tsv`
- `conservation/gene_age.tsv`
- `conservation/reports/conservation_qc.tsv`

The QC report should include:

- number of query transcripts;
- number of successfully mapped transcripts;
- number of segments before and after duplicate removal;
- number of conserved transcripts by each filter;
- number of homologous genes per species;
- distribution summaries for mapped length, match ratio, transcript coverage
  and exon/intron coverage.

## References

- LncBook 2.0 database: https://ngdc.cncb.ac.cn/lncbook/
- LncBook 2.0 article: https://pmc.ncbi.nlm.nih.gov/articles/PMC9825513/
- UCSC chain/net and pslMap-style coordinate mapping resources:
  https://genome.ucsc.edu/

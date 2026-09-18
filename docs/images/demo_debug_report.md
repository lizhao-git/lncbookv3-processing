# ClinVar variant annotation debug report

## Overview

| # | step | records | errors | warnings | status |
| --- | --- | --- | --- | --- | --- |
| 1 | `5_features_bed_report.tsv` | 7 | 1 | 0 | FAIL |
| 2 | `4_variant_bed_report.tsv` | 3214 | 0 | 1 | WARN |
| 3 | `1_annotation_validation_report.tsv` | 1471438 | 0 | 0 | OK |
| 4 | `2_clinvar_validation_report.tsv` | 8842 | 0 | 0 | OK |
| 5 | `3_clinvar_filter_report.tsv` | - | 0 | 0 | OK |
| 6 | `6_clinvar_summary.tsv` | 8 | - | - | OK |

## Errors & Warnings

### `5_features_bed_report.tsv`

- **error_count**: 1
- **errors**: line 4: BED end must be > start

### `4_variant_bed_report.tsv`

- **warning_count**: 1
- **warnings**: line 2170: BED has extra columns but fewer than BED6 columns

## Step details

### `5_features_bed_report.tsv`

| metric | value |
| --- | --- |
| records | 7 |
| chrom_count | 1 |
| error_count | 1 |
| warning_count | 0 |
| errors | line 4: BED end must be > start |

### `4_variant_bed_report.tsv`

| metric | value |
| --- | --- |
| records | 3214 |
| chrom_count | 1 |
| error_count | 0 |
| warning_count | 1 |
| warnings | line 2170: BED has extra columns but fewer than BED6 columns |

### `1_annotation_validation_report.tsv`

| metric | value |
| --- | --- |
| format | gtf |
| input_compression | plain |
| parsed_records | 1471438 |
| feature_exon | 1471438 |
| feature_gene | 214300 |
| feature_transcript | 506089 |
| feature_CDS | 5405 |
| feature_5UTR | 3120 |
| feature_3UTR | 4011 |
| error_count | 0 |
| warning_count | 0 |

### `2_clinvar_validation_report.tsv`

| metric | value |
| --- | --- |
| records | 8842 |
| chrom_count | 1 |
| error_count | 0 |
| warning_count | 0 |

### `3_clinvar_filter_report.tsv`

| metric | value |
| --- | --- |
| kept_variants | 3214 |
| label_Pathogenic | 1870 |
| label_Benign | 1042 |
| label_Risk_factor | 302 |

### `6_clinvar_summary.tsv`

| feature_type | label | count |
| --- | --- | --- |
| exon | ALL | 3214 |
| exon | Pathogenic | 1870 |
| exon | Benign | 1042 |
| intron | ALL | 214 |
| intron | Pathogenic | 93 |
| cds | ALL | 121 |
| 5UTR | ALL | 41 |
| 3UTR | ALL | 38 |


# Local Nextflow modules

Each directory under `modules/local/` is a reusable Nextflow DSL2 process with:

- `main.nf`: the process implementation.
- `meta.yml`: short nf-core-style metadata.
- `environment.yml`: conda dependencies for standalone reuse.

These modules intentionally wrap project scripts as small, composable command-line tools. Domain-level logic is in `subworkflows/local/`.

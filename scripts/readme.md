# Preprocess for CCLONE

These two scripts show how to run the variant calling with cellSNP-lite, and how to convert the output to anndata format.

## (1) Run cellSNP-lite and annotate output with vep:
run_cellSNP-lite.sh
### Requirements
- [cellsnp-lite](https://cellsnp-lite.readthedocs.io/en/latest/) installed and runable via the command line (added to $PATH - or edit path in run_cellSNP-lite.sh)
- [vep](https://www.ensembl.org/info/docs/tools/vep/index.html) installed and runable via the command line
- data

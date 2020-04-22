# Blue Lab pipeline

Tip: run R with `R -q --vanilla --args <arguments> < script.R`
instead of using `Rscript`, so R will print out a full log of all
commands and messages in your script.

Required R packages: argparser, SeqArray, SeqVarTools, SNPRelate, GENESIS

## Convert to GDS

Use the SeqArray package to convert files.

If you have plink files (BED/BIM/FAM), use `seqBED2GDS`. If you have VCF, use `seqVCF2GDS`. 

Example for VCF:
http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test_seq.html#convert-vcf-to-gds

Script for BED:
[plink_to_gds.R](plink_to_gds.R)


## QC


## Relatedness

http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test_seq.html#population-structure-and-relatedness


1. [ld_pruning.R](ld_pruning.R)

# Blue Lab pipeline

Tip: run R with `R -q --vanilla --args <arguments> < script.R`
instead of using `Rscript`, so R will print out a full log of all
commands and messages in your script.

Required R packages: argparser, SeqArray, SeqVarTools, SNPRelate, GENESIS

## Running with the PBS job scheduler

From beluga1, use [runRscript.sh](runRscript.sh). Use the `-v` flag to
pass the name of the R script (after `R=`) and any arguments to the
script (after `args=`).

```
qsub -N jobname -v R="myscript.R",args="--arg1 arg1 --arg2 arg2" runRscript.sh
```

To run a script by chromosome, use the `-t` flag to submit each
chromosome as a separate task in an array job.

```
qsub -N ld_pruning -v R="ld_pruning.R",args="--maf 0.01 --missing 0.01" -t 1-22 runRscript.sh
```


## Convert to GDS

Use the SeqArray package to convert files.

If you have plink files (BED/BIM/FAM), use `seqBED2GDS`. If you have VCF, use `seqVCF2GDS`. 

Example for VCF:
http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test_seq.html#convert-vcf-to-gds

Script for BED:
[plink_to_gds.R](plink_to_gds.R)


## QC

Outline of recommended QC steps for GWAS data:
http://bioconductor.org/packages/release/bioc/vignettes/GWASTools/inst/doc/DataCleaning.pdf

The above document uses GWASTools, but SeqVarTools can be used
instead.

### Sample QC

#### Annotated vs genetic sex check

Ideally one would use X and Y intensity values (array data) or X and Y
read depth (sequencing data) to infer genetic sex. If these are not
available, can use X chromosome heterozygosity and Y chromosome
missingness to check sex.

#### Missing call rate

[sample_qc.R](sample_qc.R)


### Variant QC

#### Missing call rate and allele frequency

[variant_qc.R](variant_qc.R)

#### Hardy-Weinberg equilibrium

When a dataset contains samples from multiple populations with
different allele frequencies, it is recommended to run HWE in each
population separately.

[hwe.R](hwe.R)


## Relatedness

http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test_seq.html#population-structure-and-relatedness


1. [ld_pruning.R](ld_pruning.R)
2. [pcs_and_grm.R](pcs_and_grm.R)

As an additional QC step, compare pairwise relatedness estimates from
PC-Relate with expected values from the pedigree. Example code:
https://github.com/UW-GAC/analysis_pipeline/blob/master/R/pedigree_check.R


## Association testing

http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test_seq.html#association-tests

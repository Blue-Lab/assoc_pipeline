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
qsub -N ld_pruning -v R="/acct/sdmorris/code/assoc_pipeline/ld_pruning.R",args="myfile.gds --maf 0.01 --missing 0.01" -t 1-22 runRscript.sh
```

To run a script in parallel, use the `-l` flag to request multiple
processors on a node, in addition to the `--num_cores` argument to the
R script.

```
qsub -N sample_qc -v R="/acct/sdmorris/code/assoc_pipeline/sample_qc.R",args="myfile.gds --num_cores 12" -l nodes=1:ppn=12 runRscript.sh
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

## Dependencies

This pipeline primarily uses R internally. It was built and tested on CentOS 7
and CentOS8 machines with R version 4.0.5, with the following package versions:

argparser\_0.7.1
Biobase\_2.50.0
BiocGenerics\_0.36.1
coxmeg_1.0.13
dplyr\_1.0.8
GENESIS\_2.20.1
gdsfmt\_1.26.1
GGally\_2.1.1
ggplot2\_3.3.5
magrittr\_2.0.2
RColorBrewer\_1.1-2
Rcpp_1.0.8.2
SeqArray\_1.30.0
SeqVarTools\_1.28.1
SNPRelate\_1.24.0
tidyr\_1.2.0

Additionally, this pipeline also has suggested dependencies PLINK v1.90b6.14
64-bit and [RUTH](https://github.com/statgen/ruth) (November 2020 release).

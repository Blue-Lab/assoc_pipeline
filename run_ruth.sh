#! /usr/bin/env bash

# Set default args
FIELD=GT
CLEANUP=false
HELP=false

while getopts ":p:d:P:o:g:D:k:v:s:n:f:N:ch" opt; do
  case ${opt} in
    p)
      # Path to phenotype ADF in .rds format
      PHENO_FILE=$OPTARG
      ;;
    d)
      # Name of DX variable in PHENO_FILE. Coded as 0/1.
      DX=$OPTARG
      ;;
    P)
      # Path to PLINK .bim/.bam/.fam files (prefix only)
      PLINK_PREF=$OPTARG
      ;;
    o)
      # Characters to prepend to output files
      OUT_PREF=$OPTARG
      ;;
    g)
      # Path to GDS file
      GDS_FILE=$OPTARG
      ;;
    D)
      # divobj for PC-AiR (.rds)
      DIV_OBJ=$OPTARG
      ;;
    k)
      # kinobj for PC-AiR (.rds) 
      KIN_OBJ=$OPTARG
      ;;
    v)
      # Vector of variant IDs to include (.rds)
      VARIANT_ID=$OPTARG
      ;;
    n)
      # Number of PCs to include fo PC-AiR
      N_PC=$OPTARG
      ;;
    f)
      # -GT option for RUTH
      FIELD=$OPTARG
      ;;
    N)
      # Number of cores to use where it is an option.
      NUM_CORE=$OPTARG
      ;;
    c)
      # Flag to delete intermediate files.
      CLEANUP=true
      ;;
    h)
      # Flag to print the following message:
      echo "\
        This script uses the bash utility getopts to parse command line
        options. For argument descriptions, check the \`while getopts\` section at 
        the head of run_ruth.sh"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

get_controls.R $PHENO_FILE $DX --out_file ${OUT_PREF}controls.rds

pcair.R $GDS_FILE $DIV_OBJ $KIN_OBJ --variant_id $VARIANT_ID \
  --sample_id ${OUT_PREF}controls.rds --num_core $NUM_CORE --out_prefix $OUT_PREF

add_fid.R $PLINK_PREF ${OUT_PREF}pcair_unrels.rds --out_file ${OUT_PREF}unrel_controls.txt

plink --recode vcf bgz \
  --bfile $PLINK_PREF --keep ${OUT_PREF}unrel_controls.txt --out ${OUT_PREF}unrel_controls

ruth_format_pcs.R ${OUT_PREF}pcair_pcs.rds ${OUT_PREF}unrel_controls.txt --out_file ${OUT_PREF}ruth_pcs.txt

ruth --vcf ${OUT_PREF}unrel_controls.vcf.gz --evec ${OUT_PREF}ruth_pcs.txt --out ${OUT_PREF}ruth.vcf \
  --field $FIELD --site-only --num-pc $N_PC

parse_ruth.R ${OUT_PREF}ruth.vcf --out_file ${OUT_PREF}ruth.rds

if [ "$CLEANUP" = true ] ; then
  rm ${OUT_PREF}controls.rds ${OUT_PREF}unrel_controls* \
    ${OUT_PREF}ruth_pcs.txt ${OUT_PREF}ruth.vcf ${OUT_PREF}pcair*
fi

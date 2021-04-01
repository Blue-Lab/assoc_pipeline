#! /usr/bin/env bash

# Set default args
FIELD=GT

while getopts ":p:d:P:o:g:D:k:v:s:n:f:N:" opt; do
  case ${opt} in
    p)
      PHENO_FILE=$OPTARG
      ;;
    d)
      DX=$OPTARG
      ;;
    P)
      PLINK_PREF=$OPTARG
      ;;
    o)
      OUT_PREF=$OPTARG
      ;;
    g)
      GDS_FILE=$OPTARG
      ;;
    D)
      DIV_OBJ=$OPTARG
      ;;
    k)
      KIN_OBJ=$OPTARG
      ;;
    v)
      VARIANT_ID=$OPTARG
      ;;
    s)
      SAMPLE_ID=$OPTARG
      ;;
    n)
      N_PC=$OPTARG
      ;;
    f)
      FIELD=$OPTARG
      ;;
    N)
      NUM_CORE=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

get_controls.R $PHENO_FILE $DX --out_prefix $OUT_PREF

pcair.R $GDS_FILE $DIV_OBJ $KIN_OBJ --variant_id $VARIANT_ID \
  --sample_id ${OUT_PREF}controls.rds --num_core $NUM_CORE --out_prefix $OUT_PREF

add_fid.R $PLINK_PREF ${OUT_PREF}pcair_unrels.rds --out_file ${OUT_PREF}unrel_controls.txt

plink --recode vcf bgz \
  --bfile $PLINK_PREF --keep ${OUT_PREF}unrel_controls.txt --out ${OUT_PREF}unrel_controls

ruth_format_pcs.R ${OUT_PREF}pcair_pcs.rds ${OUT_PREF}unrel_controls.txt --out_file ${OUT_PREF}ruth_pcs.txt

ruth --vcf ${OUT_PREF}unrel_controls.vcf.gz --evec ${OUT_PREF}ruth_pcs.txt --out ${OUT_PREF}ruth.vcf \
  --field $FIELD --site-only --num-pc $N_PC

parse_ruth.R ${OUT_PREF}ruth.vcf --out_file ${OUT_PREF}ruth.rds



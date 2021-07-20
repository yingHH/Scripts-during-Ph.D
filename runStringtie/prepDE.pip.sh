#!/usr/env bash
#
# prepDE.py "stringtie out to DESeq2 input"
#

ilst=$1
odir=$2
oname=$3

prepDE.py \
  -i $ilst \
  -g ${odir}/${oname}.gene_mtx.csv \
  -t ${odir}/${oname}.transcript_mtx.csv

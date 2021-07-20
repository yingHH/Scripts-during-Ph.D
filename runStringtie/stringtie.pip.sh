#!/usr/env bash
#
# stringtie
#

ifile=$1
gff3=$2
ogtf=$3
txt=$4

# echo $1 $2 $3 $4

stringtie \
  -e \
  -G $gff3 \
  -o $ogtf \
  -A $txt \
  $ifile

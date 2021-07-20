#!/usr/env bash
#
# samtools markdup
#

ifile=$1
odir=$2
sampName=${ifile##*/}
sampName=${sampName%%.*}
logfile=${odir}'/'${sampName}.rmdup.log


echoLog (){
  echo -e ">${sampName}: $*" | tee -a $logfile
}


rmdup (){
  ibam=$1
  odir=$2

  echoLog sort bam by name
  samtools \
    sort \
    -n \
    -o ${odir}/${sampName}.nsort.bam \
    $ibam

  echoLog mark mate tag
  samtools \
    fixmate \
    -m ${odir}/${sampName}.nsort.bam \
    ${odir}/${sampName}.fixmate.bam

  echoLog sort bam by position
  samtools \
    sort \
    -o ${odir}/${sampName}.sort.bam \
    ${odir}/${sampName}.fixmate.bam

  echoLog mark and remove duplication
  samtools \
    markdup \
    -r \
    -s \
    ${odir}/${sampName}.sort.bam \
    ${odir}/${sampName}.rmdup.bam \
    &>> $logfile

  echoLog index bam
  samtools \
    index \
    ${odir}/${sampName}.rmdup.bam

  rm \
    ${odir}/${sampName}.nsort.bam \
    ${odir}/${sampName}.fixmate.bam \
    ${odir}/${sampName}.sort.bam
}


if [ $# -eq 2 ];then
  rmdup $ifile $odir
else
  echo "Usage: rmdup.pip.sh [/path/to/read.bam] [/path/to/odir]"
fi

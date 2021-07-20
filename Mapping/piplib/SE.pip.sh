#!/bin/bash
#
# All Pipeline Pool
#
# Infile:
# Options:
# Outfile:
#
# Author: Ying Huang @HZAU, yinghuang_hzau@163.com
# Last modified: 2018-3-15
# Ver: 0.0.1

#==============================
# The following are tools Pool
#==============================

filename (){
  # $1: file to get filename
  # $2: left boundry
  # $3: right boundry
  filename=$1
  filename=${filename##*${2}}
  filename=${filename%%${3}*}
  echo $filename | tee -a ${odirName}/${sample_name}.log
}

echoLog (){
  echo -e "\n=== $* ===\n" | tee -a ${odirName}/${sample_name}.log
}

cutAdapt (){
  # $1: /path/to/read_1 $2: /path/to/read_2
  # $3: read1_name      $4: read2_name
  ifastq=$1
  ofastq=$2

  echoLog Running cutadapt ...
  ${path}cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
  -q 15,10 \
  -m 20 \
  -o $ofastq \
  $ifastq \
  2>&1 | tee -a ${odirName}/${sample_name}.log
  echoLog cutadapt Done
}

Bowtie2 (){
  # $1: /path/to/trim_read_1 $2: /path/to/trim_read_2
  # $3: out put directory for SAM file
  ifastq=$1
  osam=$2

  echoLog Running bowtie2 ...
  ${path}bowtie2 \
  -p 8 \
	-x ${genome_path} \
	-U $ifastq \
	-S $osam \
	--no-unal \
	--time \
	--omit-sec-seq \
  2>&1 | tee -a ${odirName}/${sample_name}.log
  echoLog bowtie2 Done
  rm $ifastq
}

Samtools (){
  # $1: /path/to/SAM $2: outfile.prefix
  isam=$1
  sampName=$2

  echoLog sam to bam
  ${path}samtools view \
  -b \
  -S \
  -o ${sampName}.bam \
  $isam \
  2>&1 | tee -a ${odirName}/${sample_name}.log

  echoLog sort bam by name
  samtools \
    sort \
    -n \
    -o ${sampName}.nsort.bam \
    ${sampName}.bam \
    2>&1 | tee -a ${odirName}/${sample_name}.log

  echoLog mark mate tag
  samtools \
    fixmate \
    -m ${sampName}.nsort.bam \
    ${sampName}.fixmate.bam \
    2>&1 | tee -a ${odirName}/${sample_name}.log

  echoLog sort bam by position
  samtools \
    sort \
    -o ${sampName}.sort.bam \
    ${sampName}.fixmate.bam \
    2>&1 | tee -a ${odirName}/${sample_name}.log

  echoLog mark and remove duplication
  samtools \
    markdup \
    -r \
    -s \
    ${sampName}.sort.bam \
    ${sampName}.rmdup.bam \
    2>&1 | tee -a ${odirName}/${sample_name}.log

  echoLog index bam
  samtools \
    index \
    ${sampName}.rmdup.bam \
    2>&1 | tee -a ${odirName}/${sample_name}.log

  rm \
    ${sampName}.sam \
    ${sampName}.bam \
    ${sampName}.nsort.bam \
    ${sampName}.fixmate.bam \
    ${sampName}.sort.bam
}


# pipline
SE_pipline (){
  ifastq=$1
  odir=$2

  date | tee -a ${sample_name}.log
  echoLog Using SE_pipline ...
  # do cutadapt
  cutAdapt $ifastq ${odir}/trim_${sample_name}.fastq.gz
  # do bowtie2
  Bowtie2 ${odir}/trim_${sample_name}.fastq.gz ${odir}/${sample_name}.sam
  # do samtools
  Samtools ${odir}/${sample_name}.sam ${odir}/${sample_name}
  # end time
  date | tee -a ${odir}/${sample_name}.log
  echo -e "\n******************************FIN*********************************\n" | tee -a ${odir}/${sample_name}.log
}

#======
# Main
#======

ifastqName=$1
odirName=$2
genome_path=$3
sample_name=(`filename $ifastqName / .f`)
logName=${odirName}
path=''

if [ $# -eq 3 ]
then
  SE_pipline $ifastqName $odirName
else
  echo "Usage:SE.pip.sh [/path/to/read.fastq.gz] [/path/to/odir] [/path/to/genome]"
fi

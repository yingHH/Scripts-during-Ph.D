#!/bin/bash
#
# All Pipeline Pool
#
# Infile:
# Options:
# Outfile:
#
# Author: Ying Huang @HZAU, yinghuang_hzau@163.com
# Last modified: 2017-3-27
# Ver: 0.1.2

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
  echo $filename | tee -a ${sample_name}.log
}

echoLog (){
  echo -e "\n=== $* ===\n" | tee -a ${sample_name}.log
}

fastQC (){
  # $1:[infile] $2:sample_name
  echoLog Running FastQC ...
  mkdir ./${2}_QC
  ${path}fastqc \
  --extract \
  -o ${2}_QC \
  $1 \
  2>&1 | tee -a ${sample_name}.log
  echoLog FastQC Done
}

cutAdapt (){
  # $1: /path/to/read_1 $2: /path/to/read_2
  # $3: read1_name      $4: read2_name

  # adapter base
    # ... Novogene ...
  Adapt_Novogene_R1_5='TCTTTCCCTACACGACGCTCTTCCGATCT'
  Adapt_Novogene_R2_5='AGATCGGAAGAGCGTCGTGTAGGGAAAGA'
  Adapt_Novogene_R1_3='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
  Adapt_Novogene_R2_3='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
    # ... Truseq ...
  Adapt_Truseq_R1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
  Adapt_Truseq_R2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
    # ... BioMarker ...
  Adapt_BioMarker_R1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
  Adapt_BioMarker_R2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    # ... NEB ...
  Adapt_NEB_R1_5='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'
  Adapt_NEB_R2_5='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
  Adapt_NEB_R1_3='CAAGCAGAAGACGGCATACGAGAT'
  Adapt_NEB_R2_3='ATCTCGTATGCCGTCTTCTGCTTG'


  echoLog Running cutadapt ...
  ${path}cutadapt \
  -g $Adapt_Novogene_R1_5 \
  -G $Adapt_Novogene_R2_5 \
  -a $Adapt_Novogene_R1_3 \
  -A $Adapt_Novogene_R2_3 \
  -q 15,10 \
  -m 20 \
  -o trim_${3}.fastq.gz \
  -p trim_${4}.fastq.gz \
  $1 $2 \
  2>&1 | tee -a ${sample_name}.log
  echoLog cutadapt Done
}

Bowtie2 (){
  # $1: /path/to/trim_read_1 $2: /path/to/trim_read_2
  # $3: out put directory for SAM file
  echoLog Running bowtie2 ...
  ${path}bowtie2 \
  -p 8 \
	-x ${genome_path} \
	-1 $1 \
	-2 $2 \
	-S $3 \
	--no-mixed \
	--no-discordant \
	--no-unal \
	--time \
	--omit-sec-seq \
  2>&1 | tee -a ${sample_name}.log
  echoLog bowtie2 Done
  rm $1 $2
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
  2>&1 | tee -a ${sample_name}.log

  echoLog sort bam by name
  samtools \
    sort \
    -n \
    -o ${sampName}.nsort.bam \
    ${sampName}.bam \
    2>&1 | tee -a ${sample_name}.log

  echoLog mark mate tag
  samtools \
    fixmate \
    -m ${sampName}.nsort.bam \
    ${sampName}.fixmate.bam \
    2>&1 | tee -a ${sample_name}.log

  echoLog sort bam by position
  samtools \
    sort \
    -o ${sampName}.sort.bam \
    ${sampName}.fixmate.bam \
    2>&1 | tee -a ${sample_name}.log

  echoLog mark and remove duplication
  samtools \
    markdup \
    -r \
    -s \
    ${sampName}.sort.bam \
    ${sampName}.rmdup.bam \
    2>&1 | tee -a ${sample_name}.log

  echoLog index bam
  samtools \
    index \
    ${sampName}.rmdup.bam \
    2>&1 | tee -a ${sample_name}.log

  rm \
    ${sampName}.sam \
    ${sampName}.bam \
    ${sampName}.nsort.bam \
    ${sampName}.fixmate.bam \
    ${sampName}.sort.bam
}

#======
# Main
#======

read1=$1
read2=$2
odirName=$3
genome_path=$4


read1_name=(`filename $read1 / .fastq`)
read2_name=(`filename $read2 / .fastq`)
sample_name=(`filename $read1_name / _R1`)
path=''
# C:cutAdapt, B:Bowtie2, S:Samtools

CBS_pipline (){
  date | tee -a ${sample_name}.log
  echoLog Using CBS_pipline ...
  # do cutadapt
  cutAdapt $1 $2 $read1_name $read2_name
  # do bowtie2
  Bowtie2 ./trim_$read1_name.fastq.gz ./trim_$read2_name.fastq.gz ${sample_name}.sam
  # do samtools
  Samtools ./${sample_name}.sam $sample_name
  # end time
  date | tee -a ${sample_name}.log
  echo -e "\n******************************FIN*********************************\n" | tee -a ${sample_name}.log
}

if [ $# -eq 4 ]
then
  cd $odirName
  CBS_pipline $read1 $read2
else
  echo "Usage:PE.pip.sh [read1.fastq.gz] [read2.fastq.gz] [odir] [genomePath]"
fi

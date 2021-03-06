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
  echoLog Running cutadapt ...
  ${path}cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
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
  echoLog Running samtools ...
  ${path}samtools view \
  -b \
  -S \
  -o ${2}.bam \
  $1 \
  2>&1 | tee -a ${sample_name}.log

  ${path}samtools sort \
  -o ${2}.sort.bam \
  ${2}.bam \
  2>&1 | tee -a ${sample_name}.log
  rm $1 ${2}.bam

  ${path}samtools index ${2}.sort.bam \
  2>&1 | tee -a ${sample_name}.log
  echoLog samtools Done  
}

#=================================
# The following are piplines Pool
#=================================


read1_name=(`filename $1 / .fastq`)
read2_name=(`filename $2 / .fastq`)
sample_name=(`filename $read1_name / _R1`)
path=''
genome_path='/home/yingh/genome/Gallus_5.0/bt2.index/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel'
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

if [ $# -eq 2 ]
then
  CBS_pipline $1 $2
else
  echo "Usage:myPip [/path/to/read1.fastq.gz] [/path/to/read2.fastq.gz]"
fi

#!/bin/bash
#
# All Pipeline Pool
#
# Infile:
# Options:
# Outfile:
#
# Author: Ying Huang @HZAU, yinghuang_hzau@163.com
# Last modified: 2018-4-18
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

Trimmomatic (){
  #ILLUMINACLIP: 接头和引物序列在 TruSeq3-SE 中，第一步 seed 搜索允许2个碱基错配，palindrome 比对分值阈值 30，simple clip 比对分值阈值 10，palindrome 模式允许切除的最短接头序列为 8bp（默认值），palindrome 模式去除与 R1 完全反向互补的 R2（默认去除）
  #SLIDINGWINDOW: 设置4bp大小的滑窗，最低平均质量为15
  #MINLEN: 最小read长度设置为20bp

  ifastq=$1
  ofastq=$2
  
  trim_path='/home/yingh/tools/Trimmomatic-0.38/'
  adapter='/home/yingh/tools/Trimmomatic-0.38/adapters/TruSeq2-SE.fa'

  echoLog Running trimmomatic ...
  ${trim_path}trimmomatic \
    SE \
    -threads 8 \
    $ifastq \
    $ofastq \
    ILLUMINACLIP:$adapter:2:30:10 \
    SLIDINGWINDOW:4:15 \
    MINLEN:15 \
  2>&1 | tee -a ${odirName}/${sample_name}.log
  echoLog trimmomatic Done
}

HISAT2 (){
  # $1: /path/to/trim_read_1 $2: /path/to/trim_read_2
  # $3: out put directory for SAM file
  ifastq=$1
  osam=$2

  echoLog Running hisat2 ...
  ${path}hisat2 \
  -p 8 \
	-x ${genome_path} \
	-U $ifastq \
	-S $osam \
	--no-unal \
	--time \
	--omit-sec-seq \
  2>&1 | tee -a ${odirName}/${sample_name}.log
  echoLog hisat2 Done
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
  # do Trimmomatic
  Trimmomatic $ifastq ${odir}/trim_${sample_name}.fastq.gz
  # do hisat2
  HISAT2 ${odir}/trim_${sample_name}.fastq.gz ${odir}/${sample_name}.sam
  # do samtools
  Samtools ${odir}/${sample_name}.sam ${odir}/${sample_name}
  # end time
  date | tee -a ${odir}/${sample_name}.log
  echo -e "\n******************************FIN*********************************\n" | tee -a ${odirName}/${sample_name}.log
}

#======
# Main
#======

ifastqName=$1
odirName=$2
genome_path=$3
sample_name=(`filename $ifastqName / '.f(ast)?q'`)
logName=${odirName}
path=''

if [ $# -eq 3 ]
then
  SE_pipline $ifastqName $odirName
else
  echo "Usage:SE.pip.sh [/path/to/read.fastq.gz] [/path/to/odir] [/path/to/genome]"
fi

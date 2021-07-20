# mapping
> A multi thread mapping pipline program

### Example:
Input:
```
$ python3 \
    ~/bin/mapping \ # script name
    SE_pip \        # pipline name
    4 \             # thread number
    ./bam/GC-FOXL2-RNA.info.csv \   # fastq information file
    ~/genome/Gallus_5.0/chr.bt2.index/Gallus_gallus.Gallus_gallus-5.0.dna.chr \   # bowtie2 index file
    /home/data/GC-FOXL2-RNA/ \  # input directory
    ./bam/                      # output directory
```
Output:
```
-rw-r--r--. 1 yingh poultry   29 Mar 16 13:49 BGC_0_1.log
-rw-r--r--. 1 yingh poultry 431M Mar 16 13:57 BGC_0_1.sort.bam
-rw-r--r--. 1 yingh poultry 1.2M Mar 16 13:57 BGC_0_1.sort.bam.bai
-rw-r--r--. 1 yingh poultry   29 Mar 16 13:49 BGC_0_2.log
-rw-r--r--. 1 yingh poultry 443M Mar 16 13:57 BGC_0_2.sort.bam
-rw-r--r--. 1 yingh poultry 1.2M Mar 16 13:57 BGC_0_2.sort.bam.bai
```

### Usage:
```
$ python3 ~/bin/mapping -- -h
Type:        function
String form: <function main at 0x7f3c6cd75e18>
File:        ~/bin/mapping/__main__.py
Line:        13

Usage:       mapping PIPNAME NPROC INFOCSV REFGENOME IDIR ODIR
             mapping --pipName PIPNAME --nproc NPROC --infoCsv INFOCSV --refgenome REFGENOME --idir IDIR --odir ODIR
```
### Arguements:
```
PIPNAME:  shell pipline for mapping, such as "SE_pip".
NPROC:  number of thread to use.
INFOCSV:  fastq file information in '.csv' file, used for extracting input filename.
REFGENOME:  path to genome index file.
IDIR: input directory of fastq file.
ODIR: output directory, including ".log", ".sort.bam", ".sort.bam.bai".
```

### How To Add New Pipeline To Mapping Program:
1. New pipeline must be shell script. All piplines are in 'mapping/piplib' directory. And
  the Programming Interface (PI) between "piplines" (.sh) and "mapping program" (.py) is in
  'pipABC.py', such as function 'SE_pip(runPip)'.
2. New "pipline" can only get three kind of informations from "mapping". One is "input file
  path of fastq data" ( for single-end data, there is one input file; for pair-end data,
  there are two input file. ), one is "output directory", another is "genome index file".
3. The programming Interface in 'pipABC.py' receive object 'runPip', which contained five
  attribute "pipName, ifiles, odir, refgenome, signal". These attribute can send information
  to piplines.

  ```
  graph TB
      A(__main__.py)
      B(pipCtrl.py)
      C(pipABC.py)
      D(pipline.sh)

      A --API--> B
      B --multi thread--> C
      C --"ifiles"--> D
      C --"outdir"--> D
      C --"refgenome"--> D
  ```

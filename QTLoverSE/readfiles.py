# -*- coding: utf-8 -*-
"""
Created on 2019-01-14 15:05:42
Last Modified on 2019-01-14 15:05:42

Read in files.

@Author: Ying Huang
"""
from collections import deque, Iterable
import gzip
import os
import pandas as pd
import re


#=====
# Read file methods
#=====

def read_general_file(ifile):
    """Read in generally coded file.
    
    Parameters:
        ifile: <Path> input file path.
    
    Return: <generator> each line of input file."""

    with open(ifile, 'rt') as f:
        for line in f:
            yield line

        
def read_gz(ifile):
    """Read in compressed file of gzip.
    
    Parameters:
        ifile: <Path> input file path.
    
    Return: <generator> each line of input file."""

    with gzip.open(ifile, 'rt') as f:
        for line in f:
            yield line


# collect all read file methods
def read_file(ifile):
    """Read in file in proper methods.
    
    Parameters:
        ifile: <Path> input file path.
    
    Return: <generator> each line of input file."""

    if ifile.endswith('.gz'):
        return read_gz(ifile)
    else:
        return read_general_file(ifile)


#=====
# read in gff file
#=====

# function for removing 'chr' or 'chr.' from seqname
def rm_chr(x):
    if x.lower().startswith('chr.'):
        return x[4:]
    elif x.lower().startswith('chr'):
        return x[3:]
    else:
        return x

# class for picking features from attribute
class PickInfo():

    def __init__(self, df, attr):
        self.attr = attr
        self.df = df.copy()

    def _pick_attr_by_re(self, x):
        attr_re = re.search(r'{}=([^;]+)'.format(self.attr), x)
        if attr_re:
            return attr_re.group(1)

    def pick_attr(self):
        return self.df['attribute'].map(self._pick_attr_by_re)


def read_gff3(igff, features=None):
    """Read in gff3 file.
    
    Parameters:
        ifile: <Path> input gff file path, can be '.gff, .gff3, .gff.gz'.
        features: <iter> list of features that you want to pick from gff 
            attribute column.
    
    Return: <DataFrame> dataframe of input file."""

    gff = deque()
    colname = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame', 'attribute']

    for line in read_file(igff):
        if not line.startswith('#'):
            gff.append(line.strip('\n').split('\t'))

    gff = pd.DataFrame(
        list(gff),
        columns=colname,
    )

    # remove 'chr' or 'chr.' from seqname
    gff['seqname'] = gff['seqname'].map(rm_chr)

    # pick features from attribute
    if isinstance(features, Iterable):
        for attr in features:
            gff.loc[:, attr] = PickInfo(gff, attr).pick_attr()
    
    return gff[['seqname', 'start', 'end', 'source', 'feature'] + features]


#=====
# read in VCF file
#=====
def read_vcf(ifile):
    """Read in VCF file.
    
    Parameters:
        ifile: <Path> input VCF file path, can be '.vcf, .vcf.gz'.
    
    Return: <DataFrame> dataframe of input file."""

    vcf = deque()
    colname = ['chr', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    for line in read_file(ifile):
        if not line.startswith('#'):
            vcf.append(line.strip('\n').split('\t'))

    vcf = pd.DataFrame(
        list(vcf),
        columns=colname,
    )
    
    return vcf


#
# Read SE file of ROSE
#
def read_SE(ise):

    enhs = pd.read_csv(ise, skiprows=5, sep='\t')
    se = enhs[enhs['isSuper'].astype(int) == 1]
    se = se[['REGION_ID', 'CHROM', 'START', 'STOP']]

    return se[['CHROM', 'START', 'STOP']].copy()


#
# Read ucsc_refseq
#
def read_ucsc_refseq(irefseq):

    refseq = pd.read_csv(irefseq, sep='\t')
    refseq = refseq[['name', 'chrom', 'strand', 'txStart', 'txEnd', 'name2']]
    # remove 'chr' from chromosome id
    refseq['chrom'] = refseq['chrom'].map(
        lambda x: x[3:] if x.lower().startswith('chr') else x
    )

    return refseq[['chrom', 'txStart', 'txEnd']].copy()


#
# Read chromosome size file of ucsc
#
def read_chromsize(ichrom_size):

    chromsize = pd.read_csv(ichrom_size, sep='\t')
    chromsize.columns = ['chr', 'size']
    # remove 'chr' from chromosome id
    chromsize['chr'] = chromsize['chr'].map(
        lambda x: x[3:] if x.lower().startswith('chr') else x
    )

    chromsize.loc[:, 'start'] = 0
    chromsize.loc[:, 'end'] = chromsize['size'].astype(int) - 1

    return chromsize[['chr', 'start', 'end']].copy()


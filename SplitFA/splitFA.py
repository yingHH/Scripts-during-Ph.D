# -*- coding: utf-8 -*-
"""
Created on 2018-09-12 21:24:27
Last Modified on 2018-09-12 21:24:27

Split FASTA into seperated file by Chromsome.

@Author: Ying Huang
"""
import argparse
from collections import deque
import os


def main():
    args = cmd()
    split_fa(args.ifa, args.odir)


def split_fa(ifa, odir):
    check_file(ifa)
    mkdir(odir)
    fa_gen = read_fa(ifa)
    split_fa_base(fa_gen, odir)


def cmd():

    msg = """
    Split FASTA into seperated file by Chromsome.
    """

    parser = argparse.ArgumentParser(description=msg)

    parser.add_argument(
        '-o',
        dest='odir',
        action='store',
        help="Directory to store output FASTA file."
    )

    parser.add_argument(
        '-i',
        dest='ifa',
        action='store',
        help="Input a FASTA file."
    )

    args = parser.parse_args()
    
    return args


def check_file(ipath):
    assert ipath.split('.')[-1].lower() in ('fasta', 'fa', 'fna')
    assert os.path.exists(ipath)


def mkdir(odir):
    if not os.path.exists(odir):
        os.mkdir(odir)


def read_fa(ipath):
    """Read FASTA file and return a generation."""

    with open(ipath, 'r') as f:
        for line in f:
            yield line


def split_fa_base(fa_gen, odir):
    """
    fa_gen: a generation return from "read_fa()".
    """

    tmp_chr = ''
    tmp_content = deque()

    print("> Reading FASTA file ...")
    for i, line in enumerate(fa_gen):

        if line.startswith('>'):
           
            if i != 0:
                write_fa(odir, tmp_chr, tmp_content)
            
            tmp_chr = line.strip('>').split(' ')[0]
        
        tmp_content.append(line)

    write_fa(odir, tmp_chr, tmp_content)

def write_fa(odir, tmp_chr, tmp_content):

    ofile = os.path.join(odir, tmp_chr + '.fa')
    print("> Writing file '{}'...".format(os.path.basename(ofile)))
    with open(ofile, 'w') as f:
        f.write(''.join(tmp_content))
    tmp_content.clear()           

if __name__ == '__main__':
    main()

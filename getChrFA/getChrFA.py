# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 16:50:14 2018

@author: yinghuang
"""
from collections import deque

import fire


def getChrfa(ifile=None, ofile=None):
    getChr_write(GetChr, ifile, ofile)


def getChr_write(GetChrCls, ifile, ofile):
    chr_write = GetChrCls()
    chr_write(ifile, ofile)


class GetChr:
    Store = deque()
    IsChr = False 

    def is_header_chr(self, header):
        """
        Only used for '.fa' file of Ensembl.
        """
        seqType = header.split(' ')[1]
        if seqType == 'dna:chromosome':
            self.IsChr = True
        else:
            self.IsChr = False

    def if_write(self, idx, ofile):
        if idx > 1 and self.Store:
            print("> Start writing:\n...", self.Store[0].strip('\n'))
            self.write_fa(self.Store, ofile)
            #print('> out:\n', self.Store)
            self.Store.clear()

    def write_fa(self, content, ofile):
        """
        content:deque
        """
        with open(ofile, 'at') as f:
            f.write(''.join(content))

    def __call__(self, ifile, ofile):
        print("> Start reading ...")
        with open (ifile, 'rt') as f:
            for i, line in enumerate(f, 1):
                if line[0] == '>':
                    self.if_write(i, ofile)
                    self.is_header_chr(line)

                if self.IsChr is True:
                    self.Store.append(line)

            self.if_write(i, ofile)


if __name__ == '__main__':
    fire.Fire(getChrfa)

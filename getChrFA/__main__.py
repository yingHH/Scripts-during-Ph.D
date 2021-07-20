# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 10:23:57 2018

@author: yinghuang
"""
import argparse

from getChrFA.getChrFA import getChrfa


def main():
    args = commandline()
    ifile = args.ifile
    ofile = args.ofile
    getChrfa(ifile, ofile)


def commandline():
    msg ="""
    Example:
        python3 ~/bin/splitFA \\
        -i Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \\
        -o Gallus_gallus.Gallus_gallus-5.0.dna.chr.fa
    """
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument('-i',
                        dest='ifile',
                        required=True,
                        action='store',
                        help="read in a '.fa' file.")
    parser.add_argument('-o',
                        dest='ofile',
                        required=True,
                        action='store',
                        help="output a '.fa' file.")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()

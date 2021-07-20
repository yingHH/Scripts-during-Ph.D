# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 10:47:51 2018

@author: yinghuang
"""
import argparse
import os

from pipCtrl import CtrlPip
from sumLog import sum_log


PWD = os.path.split(os.path.realpath(__file__))[0]
shPath = os.path.join(PWD, 'piplib')


def main():

    args = cmd()

    print('starting ...')
    if not os.path.exists(args.infoCsv):
        raise Exception("Path '{}' not exist.".format(args.infoCsv))
    elif not os.path.exists(args.idir):
        raise Exception("Path '{}' not exist.".format(args.idir))
    elif not os.path.exists(args.odir):
        print("Make directory '{}'.".format(args.odir))
        os.mkdir(args.odir)

    res = CtrlPip(args.pipName, args.nproc,
                  args.infoCsv, args.refgenome, args.idir, args.odir)
    print(str(res))

    sum_logFiles(idir=args.odir, ofilename='summary.log.csv')


def cmd():

    pips = lst_pips()

    msg = """
    Mapping FASTQ data in multi-threads way.
    """

    parser = argparse.ArgumentParser(description=msg)

    parser.add_argument(
        '-pn', '--pipName',
        dest='pipName',
        choices=pips,
        required=True,
        help="Choose a shell pipline for mapping."
    )

    parser.add_argument(
        '-ic', '--infoCsv',
        dest='infoCsv',
        action='store',
        required=True,
        help="Input a CSV file that contains all FASTQ files information."
    )

    parser.add_argument(
        '-rg', '--refgenome',
        dest='refgenome',
        action='store',
        required=True,
        help="Path to Bowtie2 or HISAT2 genome index file."
    )

    parser.add_argument(
        '-np', '--nproc',
        dest='nproc',
        type=int,
        default=1,
        action='store',
        required=False,
        help="Number of thread to use for running pipline. Note: in all piplines, Bowtie2 and HISAT2 use 8 threads for mapping. So the real number of threads is: 'nproc * 8'."
    )

    parser.add_argument(
        '-o', '--odir',
        dest='odir',
        action='store',
        required=True,
        help="Directory path to store all results."
    )

    parser.add_argument(
        'idir',
        action='store',
        help="Directory path for input FASTQ files."
    )

    args = parser.parse_args()

    return args


def sum_logFiles(idir=None, ofilename='summary.log.csv'):
    try:
        logDf = sum_log(idir)
        logDf.to_csv(os.path.join(idir, ofilename))
        print("> Statistic all log files and write to :", ofilename)
    except:
        print("> Failed to Statistic all log files.")

def lst_pips():
    pips = os.listdir(shPath)
    return tuple(pips)


if __name__ == '__main__':
    main()

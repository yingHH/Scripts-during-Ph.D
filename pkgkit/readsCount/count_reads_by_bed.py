# -*- coding: utf-8 -*-
"""
Created on 2018-11-26 11:13:48
Last Modified on 2018-11-26 11:13:48

Write usage of the script

@Author: Ying Huang
"""
import argparse
from collections import namedtuple
from concurrent import futures
from pybedtools import bedtool
import multiprocessing
import numpy as np
import os
import pandas as pd
import re


def separate_genome_by_chr(ibed):
    """Read in BED format file and separate it by chromsome.
    
    ibed: <
    Parameters:Path> a BED format file path.
    """



    start = 0
    end = int(chr_size)
    bin_size = int(bin_size)

    start_site_lst = np.arange(start, end, bin_size)
    end_site_lst = np.append(np.arange(start + bin_size, end, bin_size), end)

    return pd.DataFrame(
        [start_site_lst, end_site_lst],
        index=['start', 'end']
    ).T.copy()


def make_genome_bin_file(chrom_size, bin_size=1e3, if_remove_chr=False):
    """Read in chrom_size file which can be download from UCSC (eg. chrom_size of chicken: http://hgdownload.soe.ucsc.edu/goldenPath/galGal5/bigZips/galGal5.chrom.sizes), then separate chromsome by bin size and output in bed format.
    
    Parameters:
    chrom_size: <Path> a chromsome size file path.
    bin_size: <int> a bin size to separate chromsome, default is 1kb.
    if_remove_chr: <bool> remove 'chr' character from chromsome name, if choose 'True'."""

    print('<Do separate genome ...>')

    chrom_size_df = pd.read_csv(chrom_size, header=None, sep='\t')

    BinFile = namedtuple('BinFile', ['chrom', 'df'])

    for chrom, size in chrom_size_df.values:

        if if_remove_chr is True and re.match(r'chr[0-9a-z_]+', chrom, re.IGNORECASE):
            chrom = chrom[3:]

        print('> Separating chromosome: ', chrom)

        tmp_df = separate_chrom(size, bin_size)
        tmp_df['chr'] = chrom
        tmp_df['name'] = '.'
        tmp_df['score'] = '.'
        tmp_df['strand'] = '.'

        bin_file = BinFile(
            chrom,
            tmp_df[
                ['chr', 'start', 'end', 'name', 'score', 'strand']
            ].copy()
        )

        yield bin_file


def count_reads_in_each_chr(ibed, ibams, ofile):
    """Count reads in each chromsome.
    
    Parameter:
    ibed: <DataFrame> separated genome bin file.
    ibams: <list> path BAM files of reads.
    ofile: <Path> file to save result."""

    print('<Do count reads>')

    bin_bed = bedtool.BedTool.from_dataframe(ibed)
    reads_count_bed = bin_bed.multi_bam_coverage(bams=ibams)

    reads_count_df = reads_count_bed.to_dataframe()
    #print('<Output file:>', ofile)
    reads_count_df.to_csv(ofile, header=None, index=False, mode='a', sep='\t')


def multi_ctrl_bin_count(ibams, chrom_size, ofile, bin_size=1e3, if_remove_chr=False, max_process=1):
    print('<Do multi control ...>')
    # prepare a generater of bin files
    bin_files = make_genome_bin_file(chrom_size, bin_size, if_remove_chr)

    # avoid duplicated ofile path
    print('<Check ofile:>', ofile)
    init = 1
    while os.path.exists(ofile):
        ofile = '{}-{}'.format(ofile, str(init))
        init += 1
    print('<ofile path is:>', ofile)

    # Now prepared for multi-threads
    is_cancel = False
    with futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        todo = []
        chroms = []
        for chrom, bin_df in bin_files:
            future = executor.submit(
                count_reads_in_each_bin,
                bin_df,
                ibams,
                ofile
            )
            todo.append(future)
            print('Scheduled for {}: {}'.format(chrom, future))
            chroms.append(chrom)

        try:
            for chrom, future in zip(chroms, futures.as_completed(todo)):
                #print('<Count chromsome:> ', chrom)
                res = future.result()
                print('<Finished:> [{}] -> [{}]'.format(chrom, future))
        except KeyboardInterrupt:
            is_cancel = True
            for future in futures:
                future.cancel()

        if is_cancel:
            executor.shutdown()

    print('<Write output to:> ', ofile)


#==========#
#   Args   #
#==========#

def cmd():

    msg = """
    Count reads in a given bin size.
    """

    parser = argparse.ArgumentParser(description=msg)

    # Required arguements

    parser.add_argument(
        '-cs', '--chrom_size_file',
        dest='chrom_size',
        action='store',
        required=True,
        help="<Path> a chromsome size file path.(eg. chrom_size of chicken: http://hgdownload.soe.ucsc.edu/goldenPath/galGal5/bigZips/galGal5.chrom.sizes)"
    )

    parser.add_argument(
        '-of', '--ofile',
        dest='ofile',
        action='store',
        required=True,
        help="<Path> file to save result."
    )

    parser.add_argument(
        'ibams',
        nargs='+',
        action='store',
        help="<list> path BAM files of reads.(eg. 1.bam 2.bam ...)"
    )

    # Optional arguements

    parser.add_argument(
        '-bs', '--bin_size',
        dest='bin_size',
        type=int,
        default=1e3,
        action='store',
        required=False,
        help="<int> a bin size to separate chromsome, default is 1000 bp."
    )

    parser.add_argument(
        '-rc', '--remove_chr',
        dest='if_remove_chr',
        default=False,
        action='store_true',
        required=False,
        help="<bool> input 'True' or 'False', remove 'chr' character from chromsome name, if choose 'True'."
    )

    parser.add_argument(
        '-mp', '--max_process',
        dest='max_process',
        type=int,
        default=1,
        action='store',
        required=False,
        help="<int> maxmum threads used to run reads count, default is 1."
    )

    args = parser.parse_args()

    return args


#==========#
#   main   #
#==========#

def main():

    args = cmd()

    multi_ctrl_bin_count(
        args.ibams, args.chrom_size,
        args.ofile, args.bin_size,
        args.if_remove_chr, args.max_process
    )


if __name__ == "__main__":
    main()

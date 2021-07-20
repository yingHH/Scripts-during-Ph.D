# -*- coding: utf-8 -*-
"""
Created on 2019-01-15 14:15:38
Last Modified on 2019-01-15 14:15:38

Find overlaped QTLs and SE subpeaks

@Author: Ying Huang
"""
from collections import namedtuple
import pandas as pd
import numpy as np
from pybedtools.bedtool import BedTool
# self packages
from QTLoverSE.readfiles import read_gff3, read_SE
from QTLoverSE.mkOtherRegion import mk_other_regions


def qtl_in_SE(ise, irefseq, ichrom_size, igff):
    """Make DataFrame SE overlaped to QTLs.
    SE QTLs: QTLs that overlaped with SEs.
    nonSE QTLs: QTLs that overlaped with non-coding regions but not overlaped 
        with SEs.
    
    Parameters:
        ise: <Path> ROSE result file named '*_AllEnhancers.table.txt'.
        irefseq: <Path> genes annotation file (*_refseq.ucsc) in a table
            separated format.
        ichrom_size: <Path> input chrom size of ucsc.
        igff: <Path> input QTLs gff file path, can be '.gff, .gff3, .gff.gz'.
        
    Return: <DataFrame> dataframe of QTLs and SE."""

    Res = namedtuple('Res', ['qtl_ov_se', 'qtl_ov_nose'])

    qtl = read_gff3(igff, ['trait_ID', 'trait', 'Map_Type', 'FlankMarkers',
        'QTL_ID'])
    # replace "''" or 'None' to np.nan
    qtl.replace({'': np.nan, 'None': np.nan}, inplace=True)
    # only use Genome mapped QTLs
    qtl = qtl[qtl['Map_Type'] == 'Genome']
    # remove QTLs that have no positions
    qtl = qtl[~(qtl.isna()[['start', 'end']].any(axis=1))]

    qtl_bed = BedTool.from_dataframe(qtl)
    
    se_bed = BedTool.from_dataframe(
        read_SE(ise).sort_values(['CHROM', 'START', 'STOP'])
    )

    # intergenic regions defined as genome regions without SEs and genes
    intergenic_regions_bed = BedTool.from_dataframe(
        mk_other_regions(ise, irefseq, ichrom_size).sort_values(
            ['chr', 'start', 'end'])
    )

    colname = ['chr', 'start', 'end'] + qtl.columns.tolist()
    qtl_ov_se = se_bed.intersect(qtl_bed, wa=True, wb=True).to_dataframe()
    qtl_ov_se.columns = colname

    # use QTLs that not in 'qtl_ov_se' to overlaped with intergenic regions
    qtl_other_bed = BedTool.from_dataframe(
        qtl[~(qtl['QTL_ID'].astype('int').isin(set(qtl_ov_se['QTL_ID'])))])
    qtl_ov_nose = intergenic_regions_bed.intersect(
        qtl_other_bed, wa=True, wb=True).to_dataframe()
    qtl_ov_nose.columns = colname

    if set(qtl_ov_se['QTL_ID']) & set(qtl_ov_nose['QTL_ID']):
        print('<Warning:> "SE QTLs" overlaped with "nonSE QTLs", which means some "nonSE QTLs" also in "SE QTLs". For scripts, it means the QTLs that not in DataFrame "qtl_ov_se" were not removed from DataFrame "qtl".')

    return Res(qtl_ov_se, qtl_ov_nose)

    

    

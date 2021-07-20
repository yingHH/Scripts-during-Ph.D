# -*- coding: utf-8 -*-
"""
Created on 2019-01-15 16:43:39
Last Modified on 2019-01-15 16:43:39

Extract non-coding regions that are not SE or coding regions,
which called other regions here.

@Author: Ying Huang
"""
import pandas as pd
from pybedtools.bedtool import BedTool
# self packages
from QTLoverSE.readfiles import read_chromsize, read_SE, read_ucsc_refseq


def mk_other_regions(ise, irefseq, ichrom_size):
    """Extract non-coding regions that are not SE or coding regions.
    
    Parameters:
        ise: <Path> ROSE result file named '*_AllEnhancers.table.txt'.
        irefseq: <Path> genes annotation file (*_refseq.ucsc) in a table
            separated format.
        ichrom_size: <Path> input chrom size of ucsc.
            
    Return: <DataFrame> bed format DataFrame file."""

    se_bed = BedTool.from_dataframe(
        read_SE(ise).sort_values(['CHROM', 'START', 'STOP'])
    )
    refseq_bed = BedTool.from_dataframe(
        read_ucsc_refseq(irefseq).sort_values(
            ['chrom', 'txStart', 'txEnd'])
    )
    chrom_size = read_chromsize(ichrom_size)
    chrom_size_bed = BedTool.from_dataframe(
        chrom_size.sort_values([
            'chr', 'start', 'end'])
    )
    
    func_regions = chrom_size_bed.intersect([refseq_bed, se_bed]).to_dataframe()

    # get non-function regions from intergenic regions of func_region
    nonfunc_regions = pd.DataFrame()
    for chrom in set(func_regions['chrom']):
        '''
        if chrom not in set(chrom_size['chr']):
            print('<Skip chromosome: {}>'.format(chrom))
            continue
        print('<Operate chromosome: {}>'.format(chrom))
        '''
        # get the start and end site of chose chromsome
        #print(chrom_size.loc[chrom_size['chr'] == str(chrom), 'start'])
        start = pd.DataFrame(
            chrom_size.loc[chrom_size['chr'] == str(chrom), 'start'],
        )
        start.columns = ['end']
        #print(start.head())
        end = pd.DataFrame(
            chrom_size.loc[chrom_size['chr'] == str(chrom), 'end'],
        )
        end.columns = ['start']
        #print(end.head())
        tmp_df = func_regions[func_regions['chrom'] == str(chrom)].drop_duplicates(
            ['chrom', 'start', 'end']
        ).copy()
        #print(tmp_df.head())

        # use end sites column of func_regions as start sites column of 
        # nonfunc_regions, also use start sites column of func_regions as end 
        # sites column of nonfunc_regions
        tmp_nf_start_sites = pd.concat(
            [start, tmp_df[['end']]],
            axis=0,
            ignore_index=True,
        )
        
        tmp_nf_end_sites = pd.concat(
            [tmp_df[['start']], end],
            axis=0,
            ignore_index=True,
        )
        
        # concat all column of nonfunc_regions
        tmp_nf = pd.concat(
            [tmp_nf_start_sites, tmp_nf_end_sites],
            axis=1,
            ignore_index=True,
        )
        #print(tmp_nf.head())
        tmp_nf.columns = ['start', 'end']
        tmp_nf.loc[:, 'chr'] = chrom
        tmp_nf = tmp_nf[['chr', 'start', 'end']]

        # if start site of regions larger than its end site, which means it is 
        # the overlaped region of func_regions rather a intergenic region of 
        # func_regions, drop these regions.
        tmp_nf = tmp_nf[tmp_nf['start'] < tmp_nf['end']]

        # concat tmp_nf to nonfunc_regions
        if nonfunc_regions.empty:
            nonfunc_regions = tmp_nf
        else:
            nonfunc_regions = pd.concat(
                [nonfunc_regions, tmp_nf],
                axis=0,
                ignore_index=True,
            )

    return nonfunc_regions.copy()



    

    
    

    


    

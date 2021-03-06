# -*- coding: utf-8 -*-
"""
Created on 2018-12-28 11:02:35
Last Modified on 2018-12-28 11:02:35

Functions for regions operation.

@Author: Ying Huang
"""
from collections import deque
import pandas as pd


BED4_Format_Colname = ('chr', 'start', 'end', 'name')


# =====
# Merge overlaped regions
# =====

# Class contains functions to modified names.
# If a new method added to the class, name of
# the method should be added to dict 'NameFuns'. 
class NameFunc:

    Tx_symbol_ids_Lst = deque()

    @classmethod
    def simple_names_type(self, names):
        """Only for get_SE_constitute_regions(). In this function,
        regions name looks like 'SE_name:Peak_name', this method
        will return overlaped name as 'SE_name:Peak_name1-Peak_name2'.
        
        names: <iter> list of names."""

        # get SE name of overlaped regions
        se_id = names[0].split(':')[0]
        # Since FIMO required length of FASTA file header less than 512
        # characters and 'peaks_id' below will be used for FASTA file
        # header, only use the start and end peaks_id in header.
        peaks = [n.split(':')[1] for n in names]
        peaks_id = '-'.join([peaks[0], peaks[-1]])

        return '{}{}{}'.format(se_id, ':', peaks_id)

    @classmethod
    def tss_names_type(self, names):
        """Only for get_tss_constitute_regions(). In this function,
        regions name looks like 'Tx_name:Peak_name', this method
        will return overlaped name as 
        'Tx_name1-Tx_name10:Peak_name1-Peak_name10'.
        
        names: <iter> list of names."""

        # get Transcripts names of overlaped regions,
        # then create 'Tx_ids_symbol' by using first Tx_id and last Tx_id,
        # finally align 'Tx_ids_symbol' to 'Tx_ids' and store it in 
        # 'Tx_symbol_ids_Lst'.
        Tx_ids = list(set([n.split(':')[0] for n in names]))
        Tx_ids_symbol = '-'.join([Tx_ids[0], Tx_ids[-1]])
        self.Tx_symbol_ids_Lst.append(
            [Tx_ids_symbol, ','.join(Tx_ids)]
        )

        # Since FIMO required length of FASTA file header less than 512
        # characters and 'peaks_id' below will be used for FASTA file
        # header, only use the start and end peaks_id in header.
        peaks = [n.split(':')[1] for n in names]
        peaks_id = '-'.join([peaks[0], peaks[-1]])

        return '{}:{}'.format(Tx_ids_symbol, peaks_id)

    @classmethod
    def mean_value_names(self, names):
        
        """names: <iter> list of names. In which names are values."""
        
        values = [float(i) for i in names]
        mean = sum(values) / len(values)

        return str(mean)


# Dict to query all methods modified names in class 'NameFunc'
NameFuns = {
    'simple_names_type': NameFunc.simple_names_type,
    'tss_names_type': NameFunc.tss_names_type,
    'mean_value_names': NameFunc.mean_value_names,
}


# Mark overlaped and non-overlaped regions by '1' and '0'
def mark_overlaps(df, ch):
    """Mark the overlaped regions in a same number.
    Note: only can be used in a chromsome per time.
    
    df: <DataFrame> BED4 format dataframe.
    ch: <str> Chromosome symbols."""

    df.columns = BED4_Format_Colname
    df = df.loc[df['chr'] == ch].copy()

    df.sort_values(['chr', 'start', 'end'], inplace=True)

    # add last region (row) coordinate to each corresponding regions (rows)
    df['last_start'] = -1
    df.iloc[1:, 4] = df.iloc[:-1, 1].values.tolist()

    df['last_end'] = -1
    df.iloc[1:, 5] = df.iloc[:-1, 2].values.tolist()

    df['order'] = [i for i in range(len(df))]

    # Mark '0' if last and current regions are not overlaped,
    # mark '1' if last and current regions are overlaped,
    # the regions overlaption can be defined as:
    # 'current start < last start < current end ||
    # current start < last end < current end'
    df['is_overlap'] = df.apply(
        lambda x: 1 if x[1] <= x[4] <= x[2] or x[1] <= x[5] <= x[2] else 0,
        axis=1,
    )

    return df


# Merge overlaped regions by marks generated by mark_overlaps()
def merge_by_marker(df, ch, name_func=None):
    """Merge overlaped regions by DataFrame generated by mark_overlaps().
    Note: only can be used in a chromsome per time.
    
    df: <DataFrame> result of mark_overlaps().
    ch: <str> Chromosome symbols.
    name_func: <Object> a function to get necessary information from name 
        of each region. The function can only receive a iterable 
        parameter. If name_func is None, the name will be '.'."""

    df = df.loc[df['chr'] == ch].copy()

    # Merge overlaped regions to a big region.
    # - Regions overlaped their corresponding last regions are marked '1' 
    # - Regions not overlaped their corresponding last regions are marked '0',
    # If one region is marked with '0', regions followed it that marked with 
    # '1' should merged with the '0' marked region to a big region, until a
    # new '0' marked region appears.
    last_order = -1
    res = []

    overlaps_interval = df[df['is_overlap'] == 0]['order'].values.tolist() + \
        [df['order'].values.tolist()[-1] + 1]

    for order in overlaps_interval:
        if last_order == -1:
            last_order = order
            continue

        overlap_df = df[
            (df['order'] >= last_order) & (df['order'] < order)].copy()
        # Get overlaped regions information
        start = overlap_df['start'].min()
        end = overlap_df['end'].max()
        name = name_func(overlap_df['name'].values.tolist().copy()) \
            if name_func is not None else '.'
        # refresh last_order
        last_order = order

        res.append([ch, start, end, name])

    return pd.DataFrame(res, columns=BED4_Format_Colname)
            

# Main function for merging regions
def merge_overlaped_regions(df, name_func=None, Tx_symbol_ids_Lst_out_path=None):
    """Merge overlaped regions.
    
    df: <DataFrame> BED4 format data to be merged.

    Parameters: <merge_by_marker()>
    name_func: <Object> a function to get necessary information from name 
        of each region. The function can only receive a iterable 
        parameter. Can be one of ['simple_names_type', 'tss_names_type',
        'mean_value_names']."""

    assert isinstance(df, pd.DataFrame)
    assert len(df.columns) == 4

    df.columns = BED4_Format_Colname

    if name_func is not None:
        name_func = NameFuns[name_func]

    res = pd.DataFrame()
    for ch in set(df['chr'].values):
        mdf = mark_overlaps(df, ch)
        overlap_df = merge_by_marker(mdf, ch, name_func)

        if res.empty:
            res = overlap_df
        else:
            res = pd.concat([res, overlap_df], axis=0)

    # If merged TSS constitute regions, 
    # store the dict between 'Tx_ids_symbol' and 'Tx_ids' as CSV file.
    Tx_symbol_ids_Lst = NameFunc.Tx_symbol_ids_Lst
    if Tx_symbol_ids_Lst and Tx_symbol_ids_Lst_out_path:
        print('<Tx_symbol_ids_Lst>:')
        pd.DataFrame(
            list(Tx_symbol_ids_Lst), 
            columns=['Tx_ids_symbol', 'Tx_ids']
        ).to_csv(
            Tx_symbol_ids_Lst_out_path, 
            index=False,
        )

    return res.copy()


# =====
# Get TSS regions based on annotation file
# =====
def get_tss(annfile, tss_up_extent, tss_dn_extent):
    """Get TSS regions based on annotation file (UCSC tab format).
    
    Parameters:
    annfile: <DataFrame> genes annotation file in a table separated format of   
        UCSC. Result of readFiles.read_annotation().
    tss_up_extent: <int> tss up stream extent length.
    tss_dn_extent: <int> tss down stream extent length.
    """
    annfile = annfile.copy()

    tss_regions = annfile[
        ['chrom', 'txStart', 'txEnd', 'name', 'strand']]

    pos_strand_mask = tss_regions['strand'] == '+'
    neg_strand_mask = tss_regions['strand'] == '-'

    tss_regions.loc[pos_strand_mask, 'tss_up'] = \
        tss_regions.loc[pos_strand_mask, 'txStart'].astype('int').map(
        lambda x: x - tss_up_extent if x > tss_up_extent else 0
    ).astype('int')
    tss_regions.loc[pos_strand_mask, 'tss_dn'] = \
        tss_regions.loc[pos_strand_mask, 'txStart'].astype('int').map(
        lambda x: x + tss_dn_extent
    ).astype('int')

    tss_regions.loc[neg_strand_mask, 'tss_up'] = \
        tss_regions.loc[neg_strand_mask, 'txEnd'].astype('int').map(
        lambda x: x - tss_dn_extent if x > tss_dn_extent else 0
    ).astype('int')
    tss_regions.loc[neg_strand_mask, 'tss_dn'] = \
        tss_regions.loc[neg_strand_mask, 'txEnd'].astype('int').map(
        lambda x: x + tss_up_extent
    ).astype('int')
    

    res = tss_regions[['chrom', 'tss_up', 'tss_dn', 'name']].copy()
    res[['tss_up', 'tss_dn']] = res[['tss_up', 'tss_dn']].astype('int64')
    res.columns = BED4_Format_Colname

    return res.copy()


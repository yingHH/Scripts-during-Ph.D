# -*- coding: utf-8 -*-
"""
Author: Ying Huang

Date: 2021-02-02 16:07:01
Last Modified by: Ying Huang
Last Modified time: 2021-02-02 16:07:01

Description: find closest maf block according to given position.
"""
import re
import pandas as pd
from collections import deque


def idx_maf(maf):
    is_block_start = False
    res = deque()
    info = {
        'seqname': None,
        'start': 0,
        'seek_start': 0,
        'seek_end': 0,
    }
    with open(maf, 'rt') as f:
        seek_pos = 0
        for line in f:
            if line.startswith('a'):
                info['seek_start'] = seek_pos
                is_block_start = True
            elif is_block_start:
                items = re.split(' +', line)
                info['seqname'] = items[1]
                info['start'] = items[2]
                is_block_start = False
            elif line.startswith('\n'):
                info['seek_end'] = seek_pos
                res.append(info.copy())
            seek_pos += len(line)
    res = pd.DataFrame(list(res))
    res.to_csv(maf + '.idx.csv', index=False)
    return res 

def maf_close_up(maf, mafidx, ref_position, query_sp=None):
    '''mafidx: CSV file resulting from "idx_maf" '''
    mafidx_df = pd.read_csv(mafidx)
    tmp_mafidx_df = mafidx_df[mafidx_df['start'] <= ref_position]

    def get_block(df):
        cur_line = df.iloc[-1].astype(str)
        print("Read in: '{}'".format('|'.join(cur_line.values.tolist())))
        seek_pos = int(cur_line['seek_start'])
        chr_len = int(cur_line['seek_end']) - int(cur_line['seek_start'])
        with open(maf, 'rt') as f:
            f.seek(seek_pos)
            return f.read(chr_len)
    res = None
    while(True):
        if tmp_mafidx_df.empty:
            print('No closest multialignment block for "{}".'.format(query_sp))
            res = None
            break
        block = get_block(tmp_mafidx_df)
        res = block
        if query_sp is None:
            break
        else:
            if ' ' + query_sp + '.' in block:
                break
            else:
                tmp_mafidx_df = tmp_mafidx_df.iloc[:-1]
                print("{} not in current block, continue. -> ".format(query_sp), end='')

    return res

def maf_close_dn(maf, mafidx, ref_position, query_sp=None):
    '''mafidx: CSV file resulting from "idx_maf" '''
    mafidx_df = pd.read_csv(mafidx)
    tmp_mafidx_df = mafidx_df[mafidx_df['start'] >= ref_position]

    def get_block(df):
        cur_line = df.iloc[0].astype(str)
        print("Read in: '{}'".format('|'.join(cur_line.values.tolist())))
        seek_pos = int(cur_line['seek_start'])
        chr_len = int(cur_line['seek_end']) - int(cur_line['seek_start'])
        with open(maf, 'rt') as f:
            f.seek(seek_pos)
            return f.read(chr_len)
    res = None
    while(True):
        if tmp_mafidx_df.empty:
            print('No closest multialignment block for "{}".'.format(query_sp))
            res = None
            break
        block = get_block(tmp_mafidx_df)
        res = block
        if query_sp is None:
            break
        else:
            if ' ' + query_sp + '.' in block:
                break
            else:
                tmp_mafidx_df = tmp_mafidx_df.iloc[1:]
                print("{} not in current block, continue. -> ".format(query_sp), end='')

    return res
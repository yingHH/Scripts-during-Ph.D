# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 12:39:09 2021

@author: Ying
"""
from Bio import AlignIO
import os
from collections import defaultdict


def get_mseq(maf, seq_name, ofile, starts, ends):
    '''
    maf: path to MAF file
    seq_name: refrence name of species and chromosome, usually at third line of maf file.
    ofile_name: output file name, will be added '.fasta' as suffix.
    start: start site of query region.
    end: end site of query region.
    '''
    
    idx_name = maf + 'index'
    
    idx = AlignIO.MafIO.MafIndex(idx_name, maf, seq_name)
    results = idx.search(starts, ends)
    
    AlignIO.write(results, ofile, "fasta")
    
    return ofile
    
def merg_fa_rawcon_seq(ifile, ofile):
    with open(ifile, 'rt') as f:
        fa = defaultdict(list)
        
        tmp_key = ''
        for line in f:
            line = line.strip('\n')
            if line.startswith('>'):
                tmp_key = line.split('.')[0]
            else:
                fa[tmp_key].append(line)
                
        new_fa = ''.join(
            ['{}\n{}\n'.format(k, ''.join(v)) for k, v in fa.items()]
        )
        new_fa = new_fa.replace('-', '')
        
        with open(ofile, 'wt') as f:
            f.write(new_fa)
          
# ！！！需要按物种拼接maf文件，存在1个物种多个scaffold拼在一个染色体的情况
def merg_fa_align_seq(ifile, ofile, ref_sp, all_headers=None):
    if all_headers is None:
        with open(ifile, 'rt') as f:
            all_headers = set([i.split('.')[0] for i in f.read().split('\n') if i.startswith('>')])
    
    last_seq_len = 0
    last_headers = []
    with open(ifile, 'rt') as f:
        fa = defaultdict(list)
        tmp_key = None
        
        for i, line in enumerate(f):
            line = line.strip('\n')
            
            # meet header
            if line.startswith('>'):
                line = line.split('.')[0]
                tmp_key = line
                
                # meet ref header
                if line == '>' + ref_sp:
                    # not the first line
                    if i != 0:
                        # supplement missing headers seq with '-'
                        mis_headers = all_headers - set(last_headers)
                        cur_seq_len = len(''.join(fa['>' + ref_sp]))
                        if mis_headers:
                            supp_seq = '-' * (cur_seq_len - last_seq_len)
                            last_seq_len = cur_seq_len
                            [fa[h].append(supp_seq) for h in mis_headers]
                        last_seq_len = cur_seq_len
                        last_headers = []
                
                last_headers.append(line)
                        
            # meet seq
            else:
                fa[tmp_key].append(line)
                
        # meet the last line
        # supplement missing headers seq with '-'
        mis_headers = all_headers - set(last_headers)
        cur_seq_len = len(''.join(fa['>' + ref_sp]))
        if mis_headers:
            supp_seq = '-' * (cur_seq_len - last_seq_len)
            [fa[h].append(supp_seq) for h in mis_headers]
                    
        # write out
        new_fa = ''.join(
            ['{}\n{}\n'.format(k, ''.join(v)) for k, v in fa.items()]
        )
        #new_fa = new_fa.replace('-', '')
        
        with open(ofile, 'wt') as f:
            f.write(new_fa)
            
        return fa

'''
get_mseq(
    'chr1_NW_020109737v1_random.maf',
    'galGal6.chr1_NW_020109737v1_random',
    './test.new',
    1160,
    1400,
)
'''
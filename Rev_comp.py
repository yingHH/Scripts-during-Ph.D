# -*- coding: utf-8 -*-
"""
Created on 2018-10-09 19:52:21
Last Modified on 2018-10-09 19:52:21

Get reverse complement sequence

@Author: Ying Huang
"""
import fire


# get reverse complement
def rev_comp(seq):
    assert isinstance(seq, str)
    
    seq = seq.upper()
    rev_seq = seq[::-1]
    
    # count complement seq
    rev_comp_seq = rev_seq
    for base, comp_base in zip(['A', 'T', 'C', 'G'], ['t', 'a', 'g', 'c']):
        rev_comp_seq = rev_comp_seq.replace(base, comp_base)
        
    return rev_comp_seq.upper()


if __name__ == '__main__':
    fire.Fire(rev_comp)
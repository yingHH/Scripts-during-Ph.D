# -*- coding: utf-8 -*-
"""
Created on 2019-03-20 22:09:02
Last Modified on 2019-03-20 22:09:02

Functions for calculating venn subsets.

@Author: Ying Huang
"""
def count_venn2(isets):
    assert len(isets) == 2

    a, b = isets

    aINb = a & b

    aNum, bNum, abNum = \
        [len(i) for i in [a, b, aINb]]

    v11 = abNum
    v10 = aNum - abNum
    v01 = bNum - abNum

    return v10, v01, v11


def count_venn3(isets):
    assert len(isets) == 3

    a, b, c = isets
    
    v100 = len(a - b - c)
    v010 = len(b - a - c)
    v110 = len((a & b) - (a & b & c))
    
    v001 = len(c - a - b)
    v101 = len((a & c) - (a & b & c))
    v011 = len((b & c) - (a & b & c))
    v111 = len(a & b & c)

    return v100, v010, v110, v001, v101, v011, v111

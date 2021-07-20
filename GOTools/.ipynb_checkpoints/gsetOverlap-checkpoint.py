# -*- coding: utf-8 -*-
"""
Author: Ying Huang
Date: 2019-11-12 12:09:14
Last Modified by: Ying Huang
Last Modified time: 2019-11-12 12:09:14
Descriptipn: count gene sets overlap between each other
"""
import pandas as pd
from scipy.stats import hypergeom, fisher_exact
from collections import deque, namedtuple


def gset_overlap(gset, header=0, geneids_sep=';'):
    """
    Count overlaps between gene sets

    Required:
    gset:<Path> or <DataFrame>.
        If <Path>, a path to a 'csv' file. If <DataFrame>, a DataFrame object.
        File should contain 2 columns. First column is 'gene set ids', second column is 'genes in gene set'.

    Optional:
    header:<int> row number of header. If no header, input:None. Default:0
    geneids_sep:<str> separater of genes in set. Default:';'
    """

    if not isinstance(gset, pd.DataFrame):
        gset = pd.read_csv(gset, header=header)

    finalset_idx = len(gset) - 1
    res = deque()

    print('<Comparing gene sets to each other ...>')
    for curset_idx in range(finalset_idx + 1):
        if curset_idx < finalset_idx:
            for cmpset_idx in range(curset_idx + 1, finalset_idx + 1):
                Set = namedtuple('Set', ['id', 'genes'])

                curset = Set(
                    gset.iloc[curset_idx, 0],
                    set(gset.iloc[curset_idx, 1].split(geneids_sep))
                )

                cmpset = Set(
                    gset.iloc[cmpset_idx, 0],
                    set(gset.iloc[cmpset_idx, 1].split(geneids_sep))
                )

                curset_only = curset.genes - cmpset.genes
                cmpset_only = cmpset.genes - curset.genes
                overlap = curset.genes & cmpset.genes

                res.append(
                    [
                        curset.id,
                        cmpset.id,
                        ';'.join(list(curset_only)),
                        ';'.join(list(cmpset_only)),
                        ';'.join(list(overlap)),
                        len(curset_only),
                        len(cmpset_only),
                        len(overlap)
                    ]
                )

    res = pd.DataFrame(
        list(res), 
        columns=[
            'id_cmp1', 'id_cmp2', 
            'genes_cmp1_only', 'genes_cmp2_only', 'genes_overlap',
            'num_cmp1_only', 'num_cmp2_only', 'num_overlap',
        ]
    )

    def hypergeom_test(x):
        """
        x:<Serise> contained 'number of genes only in group A', 'number of genes only in group B', 'number of overlap genes'
        k: number of overlap genes
        M: number of all genes
        n: number of genes in A group
        N: number of genes peaked (aka. in B group)
        note: the result is same when numbers in 'n' and 'N' are exchanged
        """
        num_grpA = x[0]
        num_grpB = x[1]
        num_overlap = x[2]

        k = num_overlap
        M = num_grpA + num_grpB + num_overlap
        n = num_grpA + num_overlap
        N = num_grpB + num_overlap

        return hypergeom.sf(k, M, n, N)

    def fisher_test(x):
        """
        x:<Serise> contained 'number of genes only in group A', 'number of genes only in group B', 'number of overlap genes'
        note: 用右端单尾检验，判断是否在正对角线（左上->右下）富集的更多
        """
        num_grpA = x[0]
        num_grpB = x[1]
        num_overlap = x[2]

        table = [
            [num_overlap, num_grpA],
            [num_grpB, num_grpA + num_grpB]
        ]

        oddsratio, pvalue = fisher_exact(table, alternative='greater')

        return pvalue

    print('<Calculating overlap significance by fisher test ...>')
    res['pvalue'] = res[['num_cmp1_only', 'num_cmp2_only', 'num_overlap']].apply(fisher_test, axis=1)

    print('<Done>')

    return res.copy()
                

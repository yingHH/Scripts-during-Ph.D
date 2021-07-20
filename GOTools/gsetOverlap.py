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
from statsmodels.stats.multitest import multipletests

# =====
# functions
# =====
def hypergeom_test(df, total_genes):
    """
    (function) Apply hypergeom test to DataFrame.

    df: <DataFrame> 3 columns of 'number of genes only in group A', 'number of genes only in group B', 'number of overlap genes'.
    total_genes: <int> the number of input DEGs to run GO.
    
    k: number of overlap genes
    M: number of all genes
    n: number of genes in A group
    N: number of genes peaked (aka. in B group)
    note: the result is same when numbers in 'n' and 'N' are exchanged
    """
    assert len(df.columns) == 3

    df.columns = ['num_grpA', 'num_grpB', 'num_overlap']

    pval = df.apply(lambda x: 1 - hypergeom.cdf(
            x['num_overlap'], total_genes, x['num_grpA'] + x['num_overlap'], x['num_grpB'] + x['num_overlap']
        ),
        axis=1,
    )
    pval.name = 'pval'

    fdr = pd.Series(multipletests(pval, method='fdr_bh')[1], name='fdr', index=pval.index)
    sig = pd.concat([pval, fdr], axis=1)

    return sig

def fisher_test(df, total_genes):
    """
    (function) Apply fisher test to DataFrame.

    x:<Serise> contained 'number of genes only in group A', 'number of genes only in group B', 'number of overlap genes'
    note: 用右端单尾检验，判断是否在正对角线（左上->右下）富集的更多
    """
    assert len(df.columns) == 3

    df.columns = ['num_grpA', 'num_grpB', 'num_overlap']
    
    _, pval = df.apply(
        lambda x: fisher_exact(
            [
                [total_genes - x['num_grpA'] - x['num_grpB'], x['num_grpA']], 
                [x['num_grpB'], x['num_overlap']]
            ]
        )
    )
    pval.name = 'pval'

    fdr = pd.Series(multipletests(pval, method='fdr_bh')[1], name='fdr')
    sig = pd.concat([pval, fdr], axis=1)

    return sig


# =====
# Main
# =====

def gset_overlap(gset, total_genes, header=0, geneids_sep=';', test_meth='hypergeom'):
    """
    Count overlaps between gene sets

    Required:
    gset:<Path> or <DataFrame>.
        If <Path>, a path to a 'csv' file. If <DataFrame>, a DataFrame object.
        File should contain 2 columns. First column is 'gene set ids', second column is 'genes in gene set'.
    total_genes: <int> the number of input DEGs to run GO.

    Optional:
    header:<int> row number of header. If no header, input:None. Default:0
    geneids_sep:<str> separater of genes in set. Default:';'
    test_meth: <str> method use for evaluate overlap significance
    """
    test_meth_dict = {
        'hypergeom': hypergeom_test,
        'fisher': fisher_exact,
    }

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
    res['overlap / cmp1'] = res['num_overlap'] / (res['num_cmp1_only'] + res['num_overlap'])
    res['overlap / cmp2'] = res['num_overlap'] / (res['num_cmp2_only'] + res['num_overlap'])

    print('<Calculating overlap significance by {} test ...>'.format(test_meth))
    
    sig = test_meth_dict[test_meth](res[['num_cmp1_only', 'num_cmp2_only', 'num_overlap']], total_genes)
    res = pd.concat([res, sig], axis=1)

    print('<Done>')

    return res.copy()
                

# -*- coding: utf-8 -*-
"""
Created on Wed May 23 23:34:57 2018

@author: Administrator
"""
from decimal import Decimal
import fire
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylabp
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import os
import pandas as pd


Label_to_file = {'CD4 (L) VS  DP (T)': 'DP_vs_CD4',
                 'CD4 (L) VS  CD4 (T)': 'CD4SP_vs_CD4',
                 'CD4 (T) VS  DP (T)': 'DP_vs_CD4SP',
                 'CD8 (L) VS  DP (T)': 'DP_vs_CD8',
                 'CD8 (L) VS  CD8 (T)': 'CD8SP_vs_CD8',
                 'CD8 (T) VS  DP (T)': 'DP_vs_CD8SP',
                 'DP (T) VS  DN (T)': 'NP_vs_DP'}

def main(setLabels = None, method=None,  title=None, ofile=None):
    """
    main(['3D4', 'PAM', 'PK'], 'Title', ofile='./venn.png')
    """
    print(setLabels)
    print(ofile)
    setLabels = setLabels.split(',')
    ifiles = [Label_to_file[lb.strip(' ')] + '_result.csv' for lb in setLabels]
    allSet = read_set(ifiles, method)
    allSetNum = count_venn3(allSet, ofile)
    plot_venn3(allSetNum, setLabels, title, ofile)
    
    
def read_set(ifiles, method):
    f1, f2, f3 = ifiles
    g1 = pd.read_csv(f1)
    g2 = pd.read_csv(f2)
    g3 = pd.read_csv(f3)
    
    pick_method = {'deg': 'get_degs_gid',
                   'up': 'get_up_gid',
                   'dn': 'get_dn_gid'}
    
    s1 = eval(pick_method[method])(g1, f1)
    s2 = eval(pick_method[method])(g2, f2)
    s3 = eval(pick_method[method])(g3, f3)
    
    return s1, s2, s3
    

def get_degs_gid(df, f):
    gid = []
    odf = df.loc[(abs(df['log2FoldChange']) >= 2) & (df['padj'] <= 0.05)]
    f = f.rstrip('.csv') + 'degs.csv'
    odf.to_csv(f, index_label='Gene ID')
    print('> Write DEGs:', f)
    gid = odf.index.tolist()
    return set(gid)

def get_up_gid(df, f):
    gid = []
    odf = df.loc[(df['log2FoldChange'] >= 2) & (df['padj'] <= 0.05)]
    f = f.rstrip('.csv') + 'degs.up.csv'
    odf.to_csv(f, index_label='Gene ID')
    print('> Write UP:', f)
    gid = odf.index.tolist()
    return set(gid)

def get_dn_gid(df, f):
    gid = []
    odf = df.loc[(df['log2FoldChange'] <= -2) & (df['padj'] <= 0.05)]
    f = f.rstrip('.csv') + 'degs.down.csv'
    odf.to_csv(f, index_label='Gene ID')
    print('> Write DN:', f)
    gid = odf.index.tolist()
    return set(gid)
    

def count_venn3(isets, ofile):
    assert len(isets) == 3

    a, b, c = isets

    aINb, aINc, bINc = [i & j for i, j in zip([a, a, b], [b, c, c])]
    abc = a & b & c

    aNum, bNum, cNum, abNum, acNum, bcNum, abcNum = \
    [len(i) for i in [a, b, c, aINb, aINc, bINc, abc]]

    v111 = abcNum
    v110, v101, v011 = [intersect - abcNum for intersect in [abNum, acNum, bcNum]]
    v100, v010, v001 = [whole - interI - interJ + v111
                        for whole, interI, interJ in
                        zip([aNum, bNum, cNum],
                            [v110, v011, v101],
                            [v101, v110, v011])]
    
    ofile = ofile[:-3]
    pd.DataFrame([list(aINb - abc),]).T.to_csv(ofile + 'ab.csv', header=None)
    pd.DataFrame([list(aINc - abc),]).T.to_csv(ofile + 'ac.csv', header=None)
    pd.DataFrame([list(bINc - abc),]).T.to_csv(ofile + 'bc.csv', header=None)
    pd.DataFrame([list(abc),]).T.to_csv(ofile + 'abc.csv', header=None)
   
    return v100, v010, v001, v110, v101, v011, v111


def plot_venn3(allSetsNum= [], setLabels = [], title=None, ofile=None):

    assert len(allSetsNum) == 7
    assert len(setLabels) == 3
    assert isinstance(ofile, str)

    plt.figure(figsize=(6,4.5)) # set fig size

    v = venn3(subsets=(2, 2, 1, 2, 1, 1, 1), # set circle size
              set_labels=setLabels,
              alpha=0.5)

    # format label value in each patch
    vSum = sum(allSetsNum)
    v100, v010, v001, v110, v101, v011, v111 = allSetsNum
    vp100, vp010, vp001, vp110, vp101, vp011, vp111 = \
    [str(Decimal(i / vSum * 100).quantize(Decimal('0.00'))) + '%' for i in allSetsNum]

    # set color in each patch
    for label, color in zip(('100', '010', '001'), ('blue', 'green', 'yellow')):
        v.get_patch_by_id(label).set_color(color)

    # set patch text
    for patch, n, p in \
        zip(('100', '010', '001', '110', '101', '011', '111'),
            (v100, v010, v001, v110, v101, v011, v111),
            (vp100, vp010, vp001, vp110, vp101, vp011, vp111)):
        v.get_label_by_id(patch).set_text('{}\n({})'.format(str(n), p))

    # add line of circle
    c = venn3_circles(subsets=(2, 2, 1, 2, 1, 1, 1))

    plt.title(title)
    plt.savefig(ofile, dpi=300)
    plt.close()
    

if __name__=='__main__':
    fire.Fire(main)
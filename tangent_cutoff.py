# -*- coding: utf-8 -*-
"""
Created on 2018-10-29 20:01:43
Last Modified on 2018-10-29 20:01:43

Find a tangent point by Lagrange's mean value.

@Author: Ying Huang
"""
from collections import namedtuple
import pandas as pd
import numpy as np


def find_cutoff(se):
    assert isinstance(se, pd.Series)

    """se: contain one column for gene expression value, and index is Gene ID"""

    se = se.map(lambda x: float(x))
    se = se[se > 0]  # 去除表达量为0的基因，使最终得到的阈值更严格
    se = se - se.min()
    #print(se.min())

    df = pd.DataFrame(se)
    df.columns = ['gexp']
    # sort gene expression values descending
    df.sort_values(['gexp'], ascending=False, inplace=True)
    df['order'] = np.arange(len(df))

    ## 计算基因表达量在对角线上的投影
    # 将“基因表达量”作为y值，“基因的顺序”作为x值
    Y_max = df['gexp'].max()
    X_max = df['order'].max()

    # 用数据点（x,y）确定对角线，及其方向上的单元向量(ex, ey)
    dis_max = (X_max**2 + Y_max**2)**(1/2)
    ey = Y_max / dis_max
    ex = X_max / dis_max

    # 用点积计算所有数据点在对角线上的投影
    df['dio_gexp'] = df.apply(lambda x: np.dot(
        x.values, np.array([ey, ex])), axis=1)

    # 投影到对角线上的y值的最小值可认为是基因表达值分布的转折点（切点）
    cut_off_gene = df['dio_gexp'].idxmin()
    df.loc[df['order'] <= df.loc[cut_off_gene, 'order'], 'class'] = 0
    df.loc[df['order'] > df.loc[cut_off_gene, 'order'], 'class'] = 1
    cut_off = namedtuple('cut_off', ['id', 'df'])

    return cut_off(cut_off_gene, df.copy())

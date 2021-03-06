# -*- coding: utf-8 -*-
"""
Author: Ying Huang

Date: 2020-05-21 15:41:23
Last Modified by: Ying Huang
Last Modified time: 2020-05-21 15:41:23

Description: 
Find tissues special GO term
"""

import numpy as np
from scipy.stats import f_oneway
import pandas as pd
from itertools import chain
from sklearn.feature_selection import SelectFdr, f_classif
from sklearn.preprocessing import StandardScaler, RobustScaler, Normalizer, MinMaxScaler
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from collections import deque, namedtuple
import time
import pickle
import json

# self
from .tissues_spec_methods import jsd
from .panther import read_go_json

# global variable
Res = namedtuple('Res', ['go_genes_df', 'go_f_se', 'go_tukey_se'])


# sub functions
def abslog(x):
    res = 0
    if x != 0:
        sign = abs(x) / x
        res = sign * np.log2(abs(x) + 1)
    return res

def scaler_row(df, iclass):
    scaler = iclass()
    df = pd.DataFrame(scaler.fit_transform(df.T).T, index=df.index, columns=df.columns)
    return df.copy()

def preprocess_data(df):
    df_log = df.applymap(abslog)
    df_minmaxScale = scaler_row(df_log, MinMaxScaler)
    return df_minmaxScale.copy()

def save_pickle(data, file):
    with open(file, 'wb') as f:
        pickle.dump(data, f)
        
def read_pickle(file):
    with open(file, 'rb') as f:
        return pickle.load(f)

# main
def sp_gos(gene_exp, go_json, sp_tissues, anova_alpha=0.05, tukey_alpha=0.05, to_pickle=None):
    """
    Get tissues special GOs based on their genes sets tissues special value.

    Required:
    gene_exp: <DataFrame> gene ids in row names, tissues id in column names, gene expression value in cells.
    go_json: <Path> path to GO result of Panther in JSON format.
    sp_tissues: <str> tissues name of special high expressed genes sets of GOs.

    Optional:
    anova_alpha: <float> P value of threshold of anova test. Default 0.05.
    tukey_alpha: <float> P value of threshold of tukey test. Default 0.05.
    to_pickle: <Path Prefix> path prefix if write out result. Default None.
    """

    # calculate genes tissues special value
    print("<{}> calculate JSD tissues special expression value".format(time.asctime( time.localtime(time.time()))))
    gene_spv = jsd(preprocess_data(gene_exp))
    
    
    # get go terms
    go = read_go_json(go_json)

    go_genes = deque()
    tissues = gene_spv.columns

    # make GO ID -> Gene List matrix
    print("<{}> make GO genes tissues special expression matrix".format(time.asctime( time.localtime(time.time()))))
    for goid, genes in go[['ID', 'genes']].values:
        gnames = genes.split(';')
        go_genes_spv = gene_spv[gene_spv.index.isin(gnames)]

        #print(go_genes_spv)

        if go_genes_spv.empty:
            print(
                ("\t<Warning zero> : GO term '{}' did not find its genes sets expression.").format(goid)
            )
            continue

        if len(go_genes_spv) == 1:
            print("\t<Warning only one> : remove term '{}' which only have one gene expression value.".format(goid))
            continue

        mean_jsd_tissues = go_genes_spv.mean(axis=0)
        if mean_jsd_tissues[sp_tissues] != mean_jsd_tissues.max():
            continue

        rdg_num = len(gnames)
        gtg_num = len(go_genes_spv)
        if rdg_num > gtg_num:
            
            print(
                ("\t<Warning> : read in '{}' genes from '{}', get '{}' genes from gene expression matrix.").format(
                    str(rdg_num), goid, str(gtg_num)
                )
            )

        go_genes.append([
            goid, 
            go_genes_spv.T.apply(lambda x: x.tolist(), axis=1).tolist().copy(), 
            [[i] * gtg_num for i in tissues]
        ])

    go_genes_df = pd.DataFrame(list(go_genes), columns=['GO ID', 'sp_value', 'sp_grp']).set_index(['GO ID'])
    print(
        "\t<Got {} GOs from {} input GOs, after filter low expressed GO genes sets across tissues>".format(
            str(len(go_genes_df)), str(len(go))
        )
    )
    

    # f test
    print("<{}> run F-test".format(time.asctime( time.localtime(time.time()))))
    def f_test(x):
        f, p = f_oneway(*x)
        #print(*x)
        return p
        
    #print(go_genes_df.head())
    f_res = go_genes_df['sp_value'].map(f_test)
    go_f_se = f_res[f_res <= anova_alpha].copy()

    print(
        "\t<Got {} GOs from {} input GOs, after F-test P <= {}>".format(
            str(len(go_f_se)), str(len(go_genes_df)), str(anova_alpha)
        )
    )
    

    # tukey_alpha
    print("<{}> run TukeyHSD-test".format(time.asctime( time.localtime(time.time()))))
    def tukeyhsd(x):
        values = list(chain(*x['sp_value']))
        grps = list(chain(*x['sp_grp']))

        mcomp = pairwise_tukeyhsd(
            values, 
            grps,
            alpha=tukey_alpha,
        ).summary().data

        mcomp = pd.DataFrame(mcomp[1:], columns=mcomp[0])

        is_sp_high_exp = (
            mcomp[
                (mcomp[['group1', 'group2']] == sp_tissues).any(axis=1)
            ]['reject'] == True
        ).all()

        return is_sp_high_exp

    tukey_res = go_genes_df.apply(tukeyhsd, axis=1)
    go_tukey_se = tukey_res[tukey_res == True]

    print(
        "\t<Got {} GOs from {} input GOs, after TukeyHSD-test P <= {}>".format(
            str(len(go_tukey_se)), str(len(go_f_se)), str(tukey_alpha)
        )
    )

    res = Res(go_genes_df, go_f_se, go_tukey_se)

    if to_pickle is not None:
        print("<{} write to> : {}".format(time.asctime( time.localtime(time.time())), '{}.pickle'.format(to_pickle)))
        save_pickle(
            res,
            '{}.pickle'.format(to_pickle)
        )
    
    return res


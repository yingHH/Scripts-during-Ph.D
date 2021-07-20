
# -*- coding: utf-8 -*-
"""
Author: Ying Huang

Date: 2020-03-24 15:33:10
Last Modified by: Ying Huang
Last Modified time: 2020-03-24 15:33:10

Description: 

"""
import os
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
from dagofun import funcsim
from dagofun.TermSimilarity import termsim

# self
from .panther import read_go_json
from .gsetOverlap import gset_overlap


# =====
# function
# =====

def topgrp(go):
    """
    df: <DataFrame> result of read_go_json
    """
    go_topgrp = go[go.apply(lambda x: x['label'] == x['top_group'], axis=1)]
    print('<Numbers of GO Terms>: ', len(go_topgrp))

    return go_topgrp.copy()

def subtopgrp(go):
    """
    df: <DataFrame> result of read_go_json
    """
    go_subtopgrp = go[go.apply(lambda x: x['label'] == x['group'], axis=1)]
    print('<Numbers of GO Terms>: ', len(go_subtopgrp))

    return go_subtopgrp.copy()

def topgrp_goids(go):
    """
    df: <DataFrame> result of read_go_json
    """
    goids_df = go.groupby(['top_group']).apply(lambda x: x['ID'].tolist())
    return goids_df.copy()

def subgrp_goids(go):
    """
    df: <DataFrame> result of read_go_json
    """
    goids_df = go.groupby(['group']).apply(lambda x: x['ID'].tolist())
    return goids_df.copy()

def mk_overlap_mtx_pval(df, thred_fdr=0.05, thred_overcmp1=0.95, thred_overcmp2=0.95):
    res = pd.DataFrame()
    
    df = df.copy()
    #pval_scale = df['pvalue'].max() - df['pvalue'].min()
    df['fdr_rate'] = df.apply(
        lambda x: 1 if ((x['fdr'] <=thred_fdr) & ((x['overlap / cmp1'] >=thred_overcmp1) | (x['overlap / cmp2'] >=thred_overcmp2))) else -1,
        axis=1,
    )
    
    for aid in set(df['id_cmp1']) | set(df['id_cmp2']):
        cmp1ids = df[df['id_cmp1'] == aid]
        cmp2ids = df[df['id_cmp2'] == aid]
        
        if cmp1ids.empty:
            cmp1overlap = pd.Series(name=aid)
        else:
            cmp1overlap = cmp1ids['fdr_rate']
            cmp1overlap.name = aid
            cmp1overlap.rename(cmp1ids['id_cmp2'].to_dict(), inplace=True)
            
        if cmp2ids.empty:
            cmp2overlap = pd.Series(name=aid)
        else:
            cmp2overlap = cmp2ids['fdr_rate']
            cmp2overlap.name = aid
            cmp2overlap.rename(cmp2ids['id_cmp1'].to_dict(), inplace=True)
            
        cmpoverlap = pd.concat([cmp1overlap, cmp2overlap], axis=0)
        #print(cmpoverlap)
        cmpoverlap[aid] = 1
        
        res = pd.concat([res, cmpoverlap], axis=1, sort=False)
        
    return res.copy()

def mk_overlap_mtx(df, thred_fdr=0.05, overlap_rate_least=0):
    res = pd.DataFrame()

    df = df.copy()

    def limit(x):
        overlap = max(x['overlap / cmp1'], x['overlap / cmp2'])

        if (x['fdr'] <=thred_fdr) and (overlap >= overlap_rate_least):
            return overlap
        else:
            return 0

    df['overlap_rate'] = df.apply(limit, axis=1)
    
    for aid in set(df['id_cmp1']) | set(df['id_cmp2']):
        cmp1ids = df[df['id_cmp1'] == aid]
        cmp2ids = df[df['id_cmp2'] == aid]
        
        if cmp1ids.empty:
            cmp1overlap = pd.Series(name=aid)
        else:
            cmp1overlap = cmp1ids['overlap_rate']
            cmp1overlap.name = aid
            cmp1overlap.rename(cmp1ids['id_cmp2'].to_dict(), inplace=True)
            
        if cmp2ids.empty:
            cmp2overlap = pd.Series(name=aid)
        else:
            cmp2overlap = cmp2ids['overlap_rate']
            cmp2overlap.name = aid
            cmp2overlap.rename(cmp2ids['id_cmp1'].to_dict(), inplace=True)
            
        cmpoverlap = pd.concat([cmp1overlap, cmp2overlap], axis=0)
        cmpoverlap[aid] = 1
        
        res = pd.concat([res, cmpoverlap], axis=1, sort=False) # concat SE columns by columns
        res.fillna(0, inplace=True)
        #print(res)
        
    return res.copy()

def mk_goSemSim_mtx(se, func='funcsim', approach='w', method='wbma', ontology='BP', res_mtx_type='squared'):
    """
    se: if func == 'funcsim', <Series> result of {topgrp,subgrp}_goids.
        if func == 'termsim', <list> a list of GO term IDs. 
    func: <str> functions of 'dagfun' to calculate GO semantic similarity. Choose one of 'funcsim', 'termsim'.
    approach: <str> see 'approach' method of 'termsim' function in [http://web.cbio.uct.ac.za/ITGOM/adagofun/Revised_ADaGOFun_2015.pdf].
    method: <str> see 'Functional Similarity Measures' in [http://web.cbio.uct.ac.za/ITGOM/adagofun/Revised_ADaGOFun_2015.pdf].
    ontology: <str> GO type, one of 'BP', 'CC', 'MF'.
    res_mtx_type: <str> result matrix type. One of 'condensed' (combinations [math] foramt), 'squared'.
    """

    if func == 'funcsim':
        goname_ids_dict = se.to_dict()
        res = funcsim(goname_ids_dict, TargetPairs=[], ontology=ontology, measure=method, output=2)
        res = pd.DataFrame(res, columns=['goSet1', 'goSet2', method])

        if res_mtx_type == 'condensed':
            res.iloc[:, 2] = res.iloc[:, 2].map(lambda x: 0 if not isinstance(x, float) else x)
            res = res
        # convert condensed dist matix to squared dist matix
        if res_mtx_type == 'squared':
            names = list(set(res['goSet1']) | set(res['goSet2']))
            res_swop = res[['goSet2', 'goSet1', method]]
            res_swop.columns = ['goSet1', 'goSet2', method]
            res = pd.concat(
                [
                    res, 
                    # set diagonal 1
                    pd.DataFrame([names, names, [1] * len(names)], index=['goSet1', 'goSet2', method]).T,
                    res_swop,
                ],
                axis=0
            ).pivot(index='goSet1', columns='goSet2', values=method)
            # change 0 to 1 on diagonal
            res = res.applymap(lambda x: 0 if isinstance(x, str) else x)
    
    if func == 'termsim':
        assert isinstance(se, (list, tuple))
        res = termsim(se, ontology=ontology, approach=approach, output=2)
        res = pd.DataFrame(res, columns=['term1', 'term2', approach])

        if res_mtx_type == 'condensed':
            res.iloc[:, 2] = res.iloc[:, 2].map(lambda x: 0 if not isinstance(x, float) else x)
            res = res
        # convert condensed dist matix to squared dist matix
        if res_mtx_type == 'squared':
            names = list(set(res['term1']) | set(res['term2']))
            res_swop = res[['term2', 'term1', approach]]
            res_swop.columns = ['term1', 'term2', approach]
            res = pd.concat(
                [
                    res, 
                    # set diagonal 1
                    pd.DataFrame([names, names, [1] * len(names)], index=['term1', 'term2', approach]).T,
                    res_swop,
                ],
                axis=0
            ).pivot(index='term1', columns='term2', values=approach)
            # change 0 to 1 on diagonal
            res = res.applymap(lambda x: 0 if isinstance(x, str) else x)

    res = res.fillna(0)
    return res.copy()


# =====
# Main
# =====

def group_go(go_json, odir, total_genes, grp_method='topgrp', test_meth='hypergeom', thred_fdr=0.05, thred_overcmp1=0.5, thred_overcmp2=0.5, figsize=(30,15)):
    """
    Group GO result base on their gene sets overlap significance.

    Required:
    go_json: <Path> GO result from Panther in JSON format.
    odir: <Paht> output directory.
    total_genes: <int> total number of input genes for GO analysis.

    Optional:
    grp_method: <str> group GO set based on what kind tree node, must be one of 'topgrp', 'subtopgrp'. Default 'topgrp'.
    test_meth: <str> method use for evaluate overlap significance. Default 'hypergeom'.
    thred_fdr: <float> thread for overlap significance. Default 0.05.
    thred_overcmp1: <float> thread for percent of overlap to first gene set of the gene sets pair. Default 0.5.
    thred_overcmp2: <float> thread for percent of overlap to second gene set of the gene sets pair. Default 0.5.
    figsize: <tuple> figure size in width and height. Default (30, 15).
    """

    if not os.path.exists(odir):
        os.mkdir(odir)

    Grp_meth_dict = {
        'topgrp': topgrp,
        'subtopgrp': subtopgrp,
    }

    opath_go = os.path.join(odir, 'go.csv')
    opath_go_overlap = os.path.join(odir, 'go_overlap.csv')
    opath_go_overlap_pval = os.path.join(odir, 'go_overlap_pval.csv')
    opath_plot = os.path.join(odir, 'go_overlap.png')

    go = read_go_json(go_json)
    go_grp = Grp_meth_dict[grp_method](go)

    go_overlap = gset_overlap(go_grp[['label', 'genes']], total_genes, test_meth)

    go_overlap_pval = mk_overlap_mtx_pval(go_overlap, thred_fdr, thred_overcmp1, thred_overcmp2)
    lnk_array = linkage(go_overlap_pval, 'ward')


    #plt.figure(figsize=(10,15))

    fig, ax = plt.subplots(
        1, 2,
        figsize=figsize,
        #sharey='col'
    )

    dn = dendrogram(lnk_array, labels=go_overlap_pval.index, orientation='left', ax=ax[0])
    #ax0_ytks = ax[0].set_yticks(ax0_yticks)
    #ax[0].get_yaxis().set_visible(False)
    ax[0].set_yticklabels([str(i) for i in range(len(dn['ivl']))], rotation=0, fontsize='small')

    # 设置单元格边界坐标
    x, y = np.mgrid[slice(0, (len(dn['ivl']) + 1) * 10, 10),
                slice(0, (len(dn['ivl']) + 1) * 10, 10)]

    ax[1].pcolormesh(
        x, y,
        go_overlap_pval.loc[dn['ivl'], dn['ivl']],
        cmap=mpl.colors.ListedColormap(sns.color_palette(['white', '#e74c3c'])),
        edgecolors='grey', 
        linewidths=0.01,
    )

    ax0_ytks = ax[0].get_yticks()
    ax[1].get_xaxis().set_visible(False)
    ax[1].yaxis.set_ticks_position('right')
    ax1_xtks = ax[1].set_xticks(ax0_ytks)
    ax1_ytks = ax[1].set_yticks(ax0_ytks)
    ax1_ytkslb = ax[1].set_yticklabels(dn['ivl'], rotation=0, fontsize='small')

    plt.tight_layout()

    go.to_csv(opath_go, index=False)
    go_overlap.to_csv(opath_go_overlap, index=False)
    go_overlap_pval.loc[dn['ivl'], dn['ivl']].to_csv(opath_go_overlap_pval)
    plt.savefig(opath_plot, dpi=300)


def group_go_rate(go_json, odir, total_genes, grp_method='topgrp', test_meth='hypergeom', thred_fdr=0.05, thred_over=0, figsize=(30,15)):
    """
    Group GO result base on their gene sets overlap rate.

    Required:
    go_json: <Path> GO result from Panther in JSON format.
    odir: <Paht> output directory.
    total_genes: <int> total number of input genes for GO analysis.

    Optional:
    grp_method: <str> group GO set based on what kind tree node, must be one of 'topgrp', 'subtopgrp'. Default 'topgrp'.
    test_meth: <str> method use for evaluate overlap significance. Default 'hypergeom'.
    thred_fdr: <float> thread for overlap significance. Default 0.05.
    thred_over: <float> thread for max percent of overlap to gene set. Default 0.
    figsize: <tuple> figure size in width and height. Default (30, 15).
    """

    if not os.path.exists(odir):
        os.mkdir(odir)

    Grp_meth_dict = {
        'topgrp': topgrp,
        'subtopgrp': subtopgrp,
    }

    opath_go = os.path.join(odir, 'go.csv')
    opath_go_overlap = os.path.join(odir, 'go_overlap.csv')
    opath_go_overlap_rate = os.path.join(odir, 'go_overlap_rate.csv')
    opath_plot = os.path.join(odir, 'go_overlap.png')

    go = read_go_json(go_json)
    go_grp = Grp_meth_dict[grp_method](go)

    go_overlap = gset_overlap(go_grp[['label', 'genes']], total_genes, test_meth)
    go_overlap_rate = mk_overlap_mtx(go_overlap, thred_fdr, thred_over)
    lnk_array = linkage(go_overlap_rate, 'ward')


    #plt.figure(figsize=(10,15))

    fig, ax = plt.subplots(
        1, 2,
        figsize=figsize,
        #sharey='col'
    )

    dn = dendrogram(lnk_array, labels=go_overlap_rate.index, orientation='left', ax=ax[0])
    #ax0_ytks = ax[0].set_yticks(ax0_yticks)
    #ax[0].get_yaxis().set_visible(False)
    ax[0].set_yticklabels([str(i) for i in range(len(dn['ivl']))], rotation=0, fontsize='small')

    # 设置单元格边界坐标
    x, y = np.mgrid[slice(0, (len(dn['ivl']) + 1) * 10, 10),
                slice(0, (len(dn['ivl']) + 1) * 10, 10)]

    ax[1].pcolormesh(
        x, y,
        go_overlap_rate.loc[dn['ivl'], dn['ivl']],
        cmap='Blues',
        edgecolors='grey', 
        linewidths=0.01,
    )

    ax0_ytks = ax[0].get_yticks()
    ax[1].get_xaxis().set_visible(False)
    ax[1].yaxis.set_ticks_position('right')
    ax1_xtks = ax[1].set_xticks(ax0_ytks)
    ax1_ytks = ax[1].set_yticks(ax0_ytks)
    ax1_ytkslb = ax[1].set_yticklabels(dn['ivl'], rotation=0, fontsize='small')

    plt.tight_layout()

    go.to_csv(opath_go, index=False)
    go_overlap.to_csv(opath_go_overlap, index=False)
    go_overlap_rate.loc[dn['ivl'], dn['ivl']].to_csv(opath_go_overlap_rate)
    plt.savefig(opath_plot, dpi=300)


def group_go_overlap_gosim(go_json, odir, grp_method='topgrp_goids', func='funcsim', method='wbma', ontology='BP', res_mtx_type='squared', figsize=(30,15), fontsize=8):
    """
    Group GO result base on their gene sets overlap rate and GO semantic similarity.

    Required:
    go_json: <Path> GO result from Panther in JSON format.
    odir: <Paht> output directory.
    total_genes: <int> total number of input genes for GO analysis.

    Optional:
    grp_method: <str> group GO set based on what kind tree node, must be one of 'topgrp_goids', 'subgrp_goids'. Default 'topgrp_goids'.
    func: <str> functions of 'dagfun' to calculate GO semantic similarity. Choose one of 'funcsim', 'termsim'.
    method: <str> see 'Functional Similarity Measures' in [http://web.cbio.uct.ac.za/ITGOM/adagofun/Revised_ADaGOFun_2015.pdf].
    ontology: <str> GO type, one of 'BP', 'CC', 'MF'.
    figsize: <tuple> figure size in width and height. Default (30, 15).
    """

    if not os.path.exists(odir):
        os.mkdir(odir)

    Grp_meth_dict = {
        'topgrp_goids': topgrp_goids,
        'subgrp_goids': subgrp_goids,
    }

    opath_go = os.path.join(odir, 'go.csv')
    opath_gosim = os.path.join(odir, 'go_sim.csv')
    opath_plot = os.path.join(odir, 'go_sim.png')

    go = read_go_json(go_json)
    go_grp = Grp_meth_dict[grp_method](go)
    gosim_mtx = mk_goSemSim_mtx(go_grp, func=func, method=method, ontology=ontology, res_mtx_type=res_mtx_type)
    lnk_array = linkage(gosim_mtx, 'ward')


    #plt.figure(figsize=(10,15))

    fig, ax = plt.subplots(
        1, 2,
        figsize=figsize,
        #sharey='col'
    )

    dn = dendrogram(lnk_array, labels=gosim_mtx.index, orientation='left', ax=ax[0])
    #ax0_ytks = ax[0].set_yticks(ax0_yticks)
    #ax[0].get_yaxis().set_visible(False)
    ax[0].set_yticklabels([str(i) for i in range(len(dn['ivl']))], rotation=0, fontsize='small')

    # 设置单元格边界坐标
    x, y = np.mgrid[slice(0, (len(dn['ivl']) + 1) * 10, 10),
                slice(0, (len(dn['ivl']) + 1) * 10, 10)]

    ax[1].pcolormesh(
        x, y,
        gosim_mtx.loc[dn['ivl'], dn['ivl']],
        cmap='Blues',
        edgecolors='grey', 
        linewidths=0.01,
    )

    ax0_ytks = ax[0].get_yticks()
    ax[1].get_xaxis().set_visible(False)
    ax[1].yaxis.set_ticks_position('right')
    ax1_xtks = ax[1].set_xticks(ax0_ytks)
    ax1_ytks = ax[1].set_yticks(ax0_ytks)
    ax1_ytkslb = ax[1].set_yticklabels(dn['ivl'], rotation=0, fontsize=fontsize)

    plt.tight_layout()

    go.to_csv(opath_go, index=False)
    gosim_mtx.loc[dn['ivl'], dn['ivl']].to_csv(opath_gosim)
    plt.savefig(opath_plot, dpi=300)



def split_cluster(overlap_mtx, thred_dist, odir, figsize=(20,16)):
    """
    Split cluster by f_cluster.

    Required:
    overlap_mtx: <DataFrame> go group names pairs are in rows and columns, their overlap values are in each cells.
    thred_dixt: <float> thread for cutoff distance.
    odir: <Path> output directory.

    Optional:
    figsize: <tuple> figure size in width and height. Default (20, 16).
    """

    if not os.path.exists(odir):
        os.mkdir(odir)
    opath_clust = os.path.join(odir, 'clust.csv')
    opath_plot = os.path.join(odir, 'overlap_clust.heatmap.png')

    lnk_array = linkage(overlap_mtx, 'ward')


    #plt.figure(figsize=(10,15))
    fig = plt.figure(constrained_layout=False)
    spec = fig.add_gridspec(
            ncols=2, nrows=1, 
            width_ratios=[1, 4], 
            wspace=0.05,
            #left=left, right=right,
        )

    ax0 = fig.add_subplot(spec[0,0])
    ax1 = fig.add_subplot(spec[0,1])

    # plot tree
    dn = dendrogram(lnk_array, labels=overlap_mtx.index, orientation='left', ax=ax0)
    overlap_mtx.loc[:, 'clust'] = fcluster(lnk_array, thred_dist, criterion='distance')
    
    clust_df = overlap_mtx.loc[dn['ivl'], 'clust'].reset_index()
    clust_boundary = []
    for i, (n, df) in enumerate(clust_df.groupby(['clust'])):
        if i < len(clust_df.groupby(['clust'])) - 1:
            boundary = df.iloc[len(df) - 1, :].name + 1
            clust_boundary.append(boundary)

    ax0.set_yticklabels([str(i) for i in range(len(dn['ivl']))], rotation=0, fontsize='small')

    # 设置单元格边界坐标
    x, y = np.mgrid[slice(0, (len(dn['ivl']) + 1) * 10, 10),
                slice(0, (len(dn['ivl']) + 1) * 10, 10)]

    # plot heatmap
    im = ax1.pcolormesh(
        x, y,
        overlap_mtx.loc[dn['ivl'], dn['ivl']],
        cmap='Blues',
        edgecolors='grey', 
        linewidths=0.01,
    )
    # plot colorbar
    cbar_ax = inset_axes(
        ax0,
        width="4%",  # width of parent_bbox width
        height="20%",  # height of parent_bbox width
        loc='upper left',
        bbox_to_anchor=(0.01, 0.03, 1, 1),
        bbox_transform=ax0.transAxes,
        borderpad=0,
    )

    fig.colorbar(
        im,
        #orientation="horizontal",
        cax=cbar_ax,
    )

    ax0_ytks = ax0.get_yticks()
    ax1.get_xaxis().set_visible(False)
    ax1.yaxis.set_ticks_position('right')
    ax1_xtks = ax1.set_xticks(ax0_ytks)
    ax1_ytks = ax1.set_yticks(ax0_ytks)
    ax1_ytkslb = ax1.set_yticklabels(dn['ivl'], rotation=0, fontsize='small')

    sns.despine(top=True, left=True, bottom=False, right=True, ax=ax0)

    # plot tree thread
    thred_x = thred_dist
    thred_y = len(overlap_mtx) * 10
    ax0.plot((thred_x, thred_x), (0, thred_y), '--', color='red')

    # plot cluster boundary
    for i in clust_boundary:
        x = (len(overlap_mtx)) * 10
        y = i * 10

        ax1.plot((0, x), (y, y), '-', color='red')

    overlap_mtx.loc[dn['ivl'], 'clust'].to_csv(opath_clust)

    fig.set_size_inches(figsize)
    plt.savefig(opath_plot, dpi=300, bbox_inches='tight')


def split_cluster_2mtx(overlap_mtx, gosim_mtx, thred_dist, odir, figsize=(7*3,4*3), fontsize=8):
    """
    Split cluster by f_cluster.

    Required:
    {overlap_mtx, gosim_mtx}: <list> two DataFrame in a list. In each DataFrame, go group names pairs are in rows and columns, their overlap values are in each cells.
    thred_dixt: <float> thread for cutoff distance.
    odir: <Path> output directory.

    Optional:
    figsize: <tuple> figure size in width and height. Default (20, 16).
    """

    if not os.path.exists(odir):
        os.mkdir(odir)
    opath_clust = os.path.join(odir, 'clust.csv')
    opath_plot = os.path.join(odir, 'overlap_clust.heatmap.png')

    
    topgrp_mtx = pd.concat([overlap_mtx, gosim_mtx], axis=1)

    lnk_topgrp = linkage(topgrp_mtx, 'ward')
    #lnk_overlap = linkage(overlap_mtx, 'ward')
    #lnk_gosim = linkage(gosim_mtx, 'ward')


    #plt.figure(figsize=(10,15))
    fig = plt.figure(constrained_layout=False)
    fig.set_size_inches(figsize)
    spec = fig.add_gridspec(
            ncols=3, nrows=1, 
            width_ratios=[1, 3, 3], 
            #height_ratios=[1, 3],
            wspace=0.01,
            hspace=0.01,
            #left=left, right=right,
        )

    #ax00 = fig.add_subplot(spec[0,0])
    #ax01 = fig.add_subplot(spec[0,1])
    #ax02 = fig.add_subplot(spec[0,2])
    ax10 = fig.add_subplot(spec[0,0])
    ax11 = fig.add_subplot(spec[0,1])
    ax12 = fig.add_subplot(spec[0,2])

    # plot topgrp tree and group tree cluster
    dn_topgrp = dendrogram(lnk_topgrp, labels=topgrp_mtx.index, orientation='left', ax=ax10)
    topgrp_mtx.loc[:, 'clust'] = fcluster(lnk_topgrp, thred_dist, criterion='distance')
    clust_df = topgrp_mtx.loc[dn_topgrp['ivl'], 'clust'].reset_index()
    clust_boundary = []
    for i, (n, df) in enumerate(clust_df.groupby(['clust'])):
        if i < len(clust_df.groupby(['clust'])) - 1:
            boundary = df.iloc[len(df) - 1, :].name + 1
            clust_boundary.append(boundary)
    
    #ax0.set_yticklabels([str(i) for i in range(len(dn['ivl']))], rotation=0, fontsize='small')

    # plot overlap tree    
    #dn_overlap = dendrogram(lnk_overlap, labels=overlap_mtx.index, orientation='top', ax=ax01)
    
    # plot gosim tree    
    #dn_gosim = dendrogram(lnk_gosim, labels=gosim_mtx.index, orientation='top', ax=ax02)

    # 设置单元格边界坐标
    x, y = np.mgrid[slice(0, (len(dn_topgrp['ivl']) + 1) * 10, 10),
                slice(0, (len(dn_topgrp['ivl']) + 1) * 10, 10)]

    # plot overlap heatmap
    im11 = ax11.pcolormesh(
        x, y,
        overlap_mtx.loc[dn_topgrp['ivl'], dn_topgrp['ivl']],
        cmap='Blues',
        edgecolors='grey', 
        linewidths=0.01,
        vmin=0,
        vmax=1,
    )    

    # plot gosim heatmap
    im12 = ax12.pcolormesh(
        x, y,
        gosim_mtx.loc[dn_topgrp['ivl'], dn_topgrp['ivl']],
        cmap='Blues',
        edgecolors='grey', 
        linewidths=0.01,
        vmin=0,
        vmax=1,
    )

    # plot colorbar
    cbar_ax11 = inset_axes(
        ax12,
        width="20%",  # width of parent_bbox width
        height="1%",  # height of parent_bbox width
        loc='upper right',
        bbox_to_anchor=(0.24, 0.04, 1, 1),
        bbox_transform=ax12.transAxes,
        borderpad=0,
    )


    fig.colorbar(
        im11,
        orientation="horizontal",
        cax=cbar_ax11,
    )


    '''
    cbar_ax12 = inset_axes(
        ax12,
        width="20%",  # width of parent_bbox width
        height="1%",  # height of parent_bbox width
        loc='upper left',
        bbox_to_anchor=(0, 0.05, 1, 1),
        bbox_transform=ax12.transAxes,
        borderpad=0,
    )


    fig.colorbar(
        im12,
        orientation="horizontal",
        cax=cbar_ax12,
    )
    '''

    #ax00.set_xticks([])
    #ax00.set_yticks([])

    #ax01.set_xticks([])
    #ax01.set_yticks([])

    #ax02.set_xticks([])
    #ax02.set_yticks([])
    
    ax10.set_yticks([])

    ax11.set_xticks([])
    ax11.set_yticks([])

    ax12.set_xticks([])
    ax12.yaxis.set_ticks_position('right')
    ax12.set_yticks([(2 * i + 1) * 5 for i in range(len(dn_topgrp['ivl']))])
    ax12.set_yticklabels(dn_topgrp['ivl'], rotation=0, fontsize=fontsize)
    #ax0_ytks = ax0.get_yticks()
    #ax1.get_xaxis().set_visible(False)
    #ax1.yaxis.set_ticks_position('right')
    #ax1_xtks = ax1.set_xticks(ax0_ytks)
    #ax1_ytks = ax1.set_yticks(ax0_ytks)


    # set title
    ax11.set_title('Genes overlap between GOs')
    ax12.set_title('GO semantic similarity')


    sns.despine(top=True, left=True, bottom=False, right=True, ax=ax10)
    #sns.despine(top=True, left=True, bottom=True, right=True, ax=ax00)
    #sns.despine(top=True, left=True, bottom=True, right=True, ax=ax01)
    #sns.despine(top=True, left=True, bottom=True, right=True, ax=ax02)

    # plot tree thread
    thred_x = thred_dist
    thred_y = len(overlap_mtx) * 10
    ax10.plot((thred_x, thred_x), (0, thred_y), '--', color='red')

    # plot cluster boundary
    for i in clust_boundary:
        x = (len(overlap_mtx)) * 10
        y = i * 10

        ax11.plot((0, x), (y, y), '-', color='red')
        ax12.plot((0, x), (y, y), '-', color='red')

    topgrp_mtx.loc[dn_topgrp['ivl'], 'clust'].to_csv(opath_clust)

    plt.savefig(opath_plot, dpi=300, bbox_inches='tight')
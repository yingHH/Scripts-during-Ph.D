# -*- coding: utf-8 -*-
"""
Created on 2019-06-10 10:26:05
Last Modified on 2019-06-10 10:26:05

KEGG enrichment analysis by 'clusterProfiler'

@Author: Ying Huang
"""
from collections import namedtuple
import pandas as pd
import numpy as np
import rpy2.robjects as ro

# self packages
from .convert import convert_geneSymbols


def kegg_enrich(readin, has_header=False, orgkegg='hsa', orgdb='org.Hs.eg.db', cutoff=0.05, cutoff_type='pvalue', is_convert_symbol_back=True, kegg_ofile=None, convert_ofile=None):
    """Run GO enrichment analysis by 'clusterProfiler'.

    Required:
    readin: <Path> or <iter>. 
        if <Path>, input gene symbols file (in 'csv' format) path, first column 
        must store gene symbols.
        if <iter>, a list of gene symbols, each one is in 'str' format.

    Optional:
    has_header: <bool> if have header or not, default 'False'.
    orgkegg: <str> organism of KEGG. All available organisms can be find in 
        "https://www.genome.jp/kegg/catalog/org_list.html".
    orgdb: <str> genome database that used in 'clusterProfiler', all available 
        database name can be listed by '.genomeDB.lst_GOfuncR_DB' function, 
        default 'org.Hs.eg.db'.
    cutoff: <float> thread to filter desired KEGG terms. Default 0.05.
    cutoff_type: <str> thread types. Can be one of 'pvalue', 'p.adjust', 
        'qvalue'. Default 'pvalue'.
    is_convert_symbol_back: <bool> if convert gene symbols of GO term to 
        original symbols. Default 'True'.
    kegg_ofile: <Path> file path to write out KEGG enrichment result. If not 
        provide, None will be write out. Default 'None'.
    convert_ofile: <Path> file path to write out gene symbols convertation 
        result. If not provide, None will be write out. Default 'None'.
    """

    assert cutoff_type in ('pvalue', 'p.adjust', 'qvalue')

    header = None if has_header else 0
    cutoff = float(cutoff)

    # read in genes
    if isinstance(readin, str):
        # if 'readin' is 'str', it's a ifile path
        genes = pd.read_csv(readin, header=header).iloc[:, 0]
    else:
        # if 'readin' is not 'str', it's a gene symbol list
        genes = readin

    # convert gene symbols to gene ids
    symbol_ids_lst = convert_geneSymbols(genes, orgdb=orgdb)

    # convert to StrVector format
    gene_ids = ro.StrVector(symbol_ids_lst.iloc[:, 1].astype('str')).r_repr()

    # run GO enrichment by 'clusterProfiler'
    print('<Run KEGG enrichment> ...')
    kk_res = ro.r(
        """
        library(clusterProfiler)
                
        kk.res <- enrichKEGG(
            {ids}, 
            organism="{organism}", 
            pvalueCutoff=0.05, 
            pAdjustMethod="BH", 
            qvalueCutoff=0.1, 
            )
            
            data.frame(kk.res@result)
        """.format(
            ids=gene_ids,
            organism=orgkegg,
        )
    )

    # convert R DataFrame result to Py DataFrame format
    kk_res = pd.DataFrame(
        np.array(kk_res).T,
        columns=list(kk_res.colnames)
    )

    # get GO terms satisfied with cutoff
    print('<Filter GO terms by [{} < {}]> ...'.format(cutoff_type, cutoff))
    kk_res = kk_res[kk_res[cutoff_type].astype('float') < cutoff]

    # convert gene symbols of GO term to original symbols
    if is_convert_symbol_back:
        print('<Convert "ENTREZID" to "Gene Names"> ...')
        ids_symbols_dict = dict(symbol_ids_lst.iloc[:, [1, 0]].values)
        kk_res['geneID'] = kk_res['geneID'].map(
            lambda x: ';'.join(
                [ids_symbols_dict[i] for i in x.split('/')]
            )
        )

    # Write out gene symbols to gene ids list
    if convert_ofile:
        print('<Writing out to: {}>'.format(convert_ofile))
        symbol_ids_lst.to_csv(convert_ofile, index=False)
    # Write out
    if kegg_ofile:
        print('<Writing out to: {}>'.format(kegg_ofile))
        kk_res.to_csv(kegg_ofile, index=False)

    res = namedtuple('Res', ['kk_res', 'convert_res'])

    return res(kk_res, symbol_ids_lst)

# -*- coding: utf-8 -*-
"""
Created on 2019-06-08 21:42:22
Last Modified on 2019-06-08 21:42:22

GO enrichment analysis by 'clusterProfiler'

@Author: Ying Huang
"""
import os
from collections import namedtuple
import pandas as pd
import numpy as np
import rpy2.robjects as ro

# self packages
from .convert import convert_geneSymbols


def go_enrich(readin, has_header=False, orgdb='org.Hs.eg.db', go_term='BP', cutoff=0.05, cutoff_type='pvalue', is_convert_symbol=True, is_convert_symbol_back=True, go_ofile=None, convert_ofile=None):
    """Run GO enrichment analysis by 'clusterProfiler'.

    Required:
    readin: <Path> or <iter>. 
        if <Path>, input gene symbols file (in 'csv' format) path, first column 
        must store gene symbols.
        if <iter>, a list of gene symbols, each one is in 'str' format.

    Optional:
    has_header: <bool> if have header or not, default 'False'.
    orgdb: <str> genome database that used in 'clusterProfiler', all available 
        database name can be listed by '.genomeDB.lst_GOfuncR_DB' function, 
        default 'org.Hs.eg.db'.
    go_term: <str> to get GO Term you wanted, can be one of character below: 
        'BP', 'CC', 'MF'. Default 'BP'.
    cutoff: <float> thread to filter desired GO terms. Default 0.05.
    cutoff_type: <str> thread types. Can be one of 'pvalue', 'p.adjust', 
        'qvalue'. Default 'pvalue'.
    is_convert_symbol_back: <bool> if convert gene symbols of GO term to 
        original symbols. Default 'True'.
    go_ofile: <Path> file path to write out GO enrichment result. If not 
        provide, None will be write out. Default 'None'.
    convert_ofile: <Path> file path to write out gene symbols convertation 
        result. If not provide, None will be write out. Default 'None'.
    """

    assert go_term in ('BP', 'CC', 'MF')
    assert cutoff_type in ('pvalue', 'p.adjust', 'qvalue')

    header = 0 if has_header else None
    cutoff = float(cutoff)

    # read in genes
    if isinstance(readin, str):
        # if 'readin' is 'str', it's a ifile path
        genes = pd.read_csv(readin, header=header).iloc[:, 0].astype('str')
    else:
        # if 'readin' is not 'str', it's a gene symbol list
        genes = readin

    if is_convert_symbol:
        # convert gene symbols to gene ids
        symbol_ids_lst = convert_geneSymbols(genes, orgdb=orgdb)
        # convert to StrVector format
        gene_ids = ro.StrVector(
            symbol_ids_lst.iloc[:, 1].astype('str')).r_repr()
    else:
        gene_ids = ro.StrVector(genes).r_repr()

    # run GO enrichment by 'clusterProfiler'
    print('<Run GO enrichment> ...')
    go_res = ro.r(
        """
        library(clusterProfiler)
        library({orgdb})
        
        go.res = enrichGO(
            gene = {igenes},
            OrgDb = {orgdb},
            keyType = 'ENTREZID',
            ont  = "{goTerm}",
            pAdjustMethod = "BH",
            pvalueCutoff  = 0.01,
            qvalueCutoff  = 0.05,
        )
        
        save(go.res, file='{save_go_res}') 
        data.frame(go.res@result)
        """.format(
            igenes = gene_ids,
            orgdb=orgdb,
            goTerm=go_term,
            save_go_res='{}.Rdata'.format(go_ofile),
        )
    )

    # convert R DataFrame result to Py DataFrame format
    go_res = pd.DataFrame(
        np.array(go_res).T,
        columns=list(go_res.colnames)
    )

    # get GO terms satisfied with cutoff
    print('<Filter GO terms by [{} < {}]> ...'.format(cutoff_type, cutoff))
    go_res = go_res[go_res[cutoff_type].astype('float') < cutoff]

    # convert gene symbols of GO term to original symbols
    if is_convert_symbol_back:
        print('<Convert "ENTREZID" to "Gene Names"> ...')
        ids_symbols_dict = dict(symbol_ids_lst.iloc[:, [1, 0]].values)
        go_res['geneID'] = go_res['geneID'].map(
            lambda x: ';'.join(
                [ids_symbols_dict[i] for i in x.split('/')]
            )
        )

    # Write out gene symbols to gene ids list
    if convert_ofile:
        print('<Writing out to: {}>'.format(convert_ofile))
        symbol_ids_lst.to_csv(convert_ofile, index=False)
    # Write out
    if go_ofile:
        print('<Writing out to: {}>'.format(go_ofile))
        go_res.to_csv(go_ofile, index=False)

    res = namedtuple('Res', ['go_res', 'convert_res'])
    
    if is_convert_symbol:
        return res(go_res, symbol_ids_lst)
    else:
        return go_res
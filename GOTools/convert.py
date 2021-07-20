# -*- coding: utf-8 -*-
"""
Created on 2019-06-08 22:01:52
Last Modified on 2019-06-08 22:01:52

Functions for converting

@Author: Ying Huang
"""
import pandas as pd
import numpy as np
import rpy2.robjects as ro


def convert_geneSymbols(symbol, fromType="SYMBOL", toType="ENTREZID", orgdb='org.Hs.eg.db'):
    """Convert gene symbols to other format('SYMBOL', 'ENTREZID' ...). Symbols 
    format can be find by `keytype(orgdb)`.

    Required:
    symbol: <iter> a list of gene symbols, each one is in 'str' format.

    Optional:
    fromType: <str> inputted gene symbols format name.
    toType: <str> wanted gene symbols format name.
    orgdb: <str> name of gene annotation database for a special species of 
        bioconducter (http://bioconductor.org/packages/3.9/BiocViews.html#___Organism).
    """

    res = ro.r(
        """
        library(clusterProfiler)
        library({orgdb})
        
        res = bitr(
            {symbol}, fromType="{fromType}", toType="{toType}", OrgDb="{orgdb}")
        """.format(
            symbol=ro.StrVector(symbol).r_repr(),
            fromType=fromType,
            toType=toType,
            orgdb=orgdb
        )
    )

    # convert R DataFrame to Py DataFrame
    res = pd.DataFrame(np.array(res).T, columns=res.colnames)

    return res.copy()


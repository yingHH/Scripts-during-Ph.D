# -*- coding: utf-8 -*-
"""
Created on 2019-06-06 15:20:01
Last Modified on 2019-06-06 15:20:01

Annotate a set of genes by GO Term using 'GOfuncR'

@Author: Ying Huang
"""
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr


def go_ann(ifile, rm_header=False, db='Homo.sapiens', pick_go_term='all', ofile=None):
    """ Map genes GO term.
    iflie: <Path> input gene symbols file (in 'csv' format) path, first column 
    must store gene symbols.
    rm_header: <bool> is remove header or not, default 'False'.
    db: <str> genome database that used in 'GOfuncR', all available database
        name can be listed by '.genomeDB.lst_GOfuncR_DB' function, default 
        'Homo.sapiens'.
    pick_go_term: <str> to filter GO Term you wanted, can be one of 'all',
        'cellular_component', 'biological_process', 'molecular_function'.
        Default 'all'.
    ofile: <Path> output file path. If not provide, None will be write 
        out. Default 'None'.
    """
    assert pick_go_term in (
        'all', 'cellular_component', 'biological_process', 'molecular_function')

    header = None if rm_header else 0

    # read in genes
    genes = pd.read_csv(ifile, header=header)
    # convert to StrVector format
    genes = robjects.StrVector(genes.iloc[:, 0])

    # annotate genes
    # import 'GOfuncR'
    gofunc = importr('GOfuncR')
    # run 'get_anno_categories'
    ann_res = gofunc.get_anno_categories(genes, database=db)

    # convert result to DataFrame format
    res = pd.DataFrame(
        np.array(ann_res).T, 
        columns=list(ann_res.colnames)
        )

    # pick go term
    if pick_go_term != 'all':
        res = res[res['domain'] == pick_go_term]
    
    # Write out
    if ofile:
        res.to_csv(ofile, index=False)

    return res.copy()

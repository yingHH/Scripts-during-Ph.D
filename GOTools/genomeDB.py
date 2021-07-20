# -*- coding: utf-8 -*-
"""
Created on 2019-06-06 17:33:01
Last Modified on 2019-06-06 17:33:01

Configure genome database

@Author: Ying Huang
"""
import pandas as pd


# =====
# Configure genome database of 'GOfuncR'
# =====

installed_DB = {
    'hg': 'Homo.sapiens',
    'mus': 'Mus.musculus',
    'hg.Org': 'org.Hs.eg.db',
    'mus.Org': 'org.Mm.eg.db',
    'gal': 'org.Gg.eg.db',
    'bv': 'org.Bt.eg.db',
}

def lst_GOfuncR_DB(is_show=True):
    """List avialable database of 'GOfuncR'. 
    """
    db_lst = pd.DataFrame(
        list(installed_DB.items()),
        columns=['species', 'database']
        )

    db_lst.set_index(['species'], drop=True, inplace=True)

    if is_show:
        print(db_lst)
    
    return db_lst.copy()

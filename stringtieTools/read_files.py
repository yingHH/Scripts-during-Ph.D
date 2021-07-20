# -*- coding: utf-8 -*-
"""
Created on 2018-12-20 17:06:15
Last Modified on 2018-12-20 17:06:15

Functions to read in files

@Author: Ying Huang
"""
from collections import namedtuple
import pandas as pd
import re


def read_transcripts_exp(ifile):
    """Read in expression file (GTF format) writed by StringTie.
    
    ifile: <Path> path to a transcripts expression file in GTF format."""

    idf = pd.read_csv(ifile, sep='\t', skiprows=2,
                      low_memory=False, header=None)

    idf = idf[(idf[1] == 'StringTie') & (idf[2] == 'transcript')]

    #features = idf[8].str.split(';', expand=True)
    #features.drop([6], axis=1, inplace=True)
    #print(features.head())

    #features_nan = features[features.isna().any(axis=1)]
    #if len(features_nan) > 0:
    #    print('<Rows contain NaN:>\n', features_nan)

    def get_values(x):

        def re_search(regex, item):
            res = re.search(regex, item)
            try:
                if res is None:
                    return '.'
                else:
                    return res.group(1)
            except:
                raise Exception('<item: {}, re_res: {}>'.format(item, res))

        '''res = namedtuple('Res', [
            'gene_id', 'transcript_id', 'ref_gene_name',
            'cov', 'FPKM', 'TPM',
        ])'''

        return ';'.join(
            list(
                re_search(r'gene_id "([^\"]+)"', x),
                re_search(r'transcript_id "([^\"]+)"', x),
                re_search(r'ref_gene_name "([^\"]+)"', x),
                re_search(r'cov "([^\"]+)"', x),
                re_search(r'FPKM "([^\"]+)"', x),
                re_search(r'TPM "([^\"]+)"', x),
            )
        )

    features = idf[8].map(
        get_values, 
    ).str.split(
        ';', expand=True
    )

    res = features

    res.columns = [
        'gene_id', 'transcript_id', 'ref_gene_name', 
        'COV', 'FPKM', 'TPM',
    ]

    return res.copy()


def read_genes_exp(ifile):
    """Read in expression file (TXT format) writed by StringTie.
    
    ifile: <Path> path to a genes expression file in TXT format."""

    return pd.read_table(ifile)[
        ['Gene ID', 'Gene Name', 'Coverage', 'FPKM', 'TPM']
    ].copy()

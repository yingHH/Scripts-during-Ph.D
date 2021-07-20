# -*- coding: utf-8 -*-
"""
Created on 2018-12-21 21:16:11
Last Modified on 2018-12-21 21:16:11

Pick desired information from GFF file, and convert it to DataFrame format.

@Author: Ying Huang
"""
from collections import defaultdict
import pandas as pd
import re


#======
#   Annotate by Refseq GFF3 file of NCBI
#======

def ncbi_rna_id_lst(igff):
    """Make list of RNA refseq ids and accession ids.
    
    igff: <Path> path to NCBI RefSeq file in GFF foramt."""

    rna_id_dict = defaultdict(list)

    with open(igff, 'rt') as f:
        for line in f:
            rna = re.search(r'ID=(rna[0-9]+);', line)
            ids = re.search(r'Genbank:([A-Z]{2}_[0-9]+\.[0-9]+)', line)
            if rna and ids:
                rna_id_dict[rna.group(1)].append(ids.group(1))

    rna_id_df = pd.DataFrame(rna_id_dict).T.reset_index()
    rna_id_df.columns = ['mask_id', 'general_id']

    return rna_id_df.copy()


def ncbi_gene_id_lst(igff):
    """Make list of gene refseq ids and Gene ids.
    
    igff: <Path> path to NCBI RefSeq file in GFF foramt."""

    gene_id_dict = defaultdict(list)

    with open(igff, 'rt') as f:
        for line in f:
            gene = re.search(r'ID=(gene[0-9]+);', line)
            ids = re.search(r'GeneID:([0-9]+)', line)
            if gene and ids:
                gene_id_dict[gene.group(1)].append(ids.group(1))

    gene_id_df = pd.DataFrame(gene_id_dict).T.reset_index()
    gene_id_df.columns = ['mask_id', 'general_id']

    return gene_id_df.copy()


def ncbi_geneid_name_generalid_lst(igff):
    """Make list of gene refseq ids and Gene names, Gene ids.
    
    igff: <Path> path to NCBI RefSeq file in GFF foramt."""

    gene_id_dict = defaultdict(list)

    with open(igff, 'rt') as f:
        for line in f:
            gene = re.search(r'ID=(gene[^;]+);', line)
            names = re.search(r'Name=([^;]+);', line)
            ids = re.search(r'GeneID:([0-9]+)', line)
            if gene and ids and names:
                gene_id_dict[gene.group(
                    1)] += [ids.group(1), names.group(1)]
    # make DataFrame, columns are 'mask_id', 'gene_name', 'general_id'
    gene_id_df = pd.DataFrame(gene_id_dict).T.reset_index()
    gene_id_df.columns = ['mask_id', 'general_id', 'Gene Name']

    return gene_id_df.copy()


def ncbi_add_generalID_to_maskID(id_type, idf, id_colname, igff):
    """Add GeneBank IDs ("NM_428001.5, XM_015284588.1, ...") to DataFrame which
    contains rna IDs ("rna0, rna1, ...").
    
    id_type: <str> only 'gene', 'tx', 'gene_name_ids' can be choice, which are      related to  function 'ncbi_gene_id_lst', 'ncbi_rna_id_lst',                 'ncbi_geneid_name_generalid_lst'.
    idf: <DataFrame> a DataFrame containing rna IDs that GeneBank IDs will be 
        add to.
    id_colname: <str | int> name of rna ID column. If 'idf' has header, the     name is a <str>, else it is a <int>.
    igff: <Path> path to NCBI RefSeq file in GFF foramt."""

    assert isinstance(idf, pd.DataFrame)
    
    if id_colname not in idf.columns.tolist():
        raise Exception(
            "'{}' is not in [{}].".format(
                id_colname, ', '.join(idf.columns.tolist())
            ) 
        )

    print('<Input DataFrame is:>\n', idf.head())
    
    print('<Making list of RNA IDs and GeneBank IDs>')\

    if id_type == 'tx':
        id_lst = ncbi_rna_id_lst(igff)
    elif id_type == 'gene':
        id_lst = ncbi_gene_id_lst(igff)
    elif id_type == 'gene_name_ids':
        id_lst = ncbi_geneid_name_generalid_lst(igff)

    print('<Merge DataFrames by column "{}">'.format(id_colname))

    if id_colname != 'mask_id':
        id_lst.rename({'mask_id': id_colname}, axis=1, inplace=True)
        
    res = pd.merge(
        id_lst, idf,
        on=id_colname,
        how='inner'
    )

    print(
        '<Input {} rows and output {} rows, since some Gene refseq IDs have no Gene IDs.>'.format(
            len(idf), len(res)
        )
    )

    return res.copy()

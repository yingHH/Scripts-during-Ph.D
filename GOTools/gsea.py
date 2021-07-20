# -*- coding: utf-8 -*-
"""
Created on 2019-09-05 14:52:40
Last Modified on 2019-09-05 14:52:40

Packages for GSEA by using ClusterProfile

@Author: Ying Huang
"""
from collections import namedtuple
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from tempfile import NamedTemporaryFile
import os
from PIL import Image

# self packages
from .convert import convert_geneSymbols


# methods to calculate values for ranking
class Rnk_Score:
    """methods to calculate values for ranking.
    Note: all methods must have three arguments 'df, idx_trt, idx_con'."""

    @classmethod
    def genelist(self, df, idx_trt, idx_con):
        return df.copy()

    @classmethod
    def Signal2Noise(self, df, idx_trt, idx_con):
        """Count the difference of means scaled by the standard deviation.

        df: <DataFrame> index column must store gene symbols, the others 
            should be numerical values(eg. gene expression values).
        idx_trt: <list> a list of column index separated of data in treatment 
            group.
        idx_con: <list> a list of column index separated of data in control
            group.
        """
        def Signal2Noise_method(se_trt, se_con):
            if se_trt.std() + se_con.std() == 0:
                res = 0
            else:
                res = (se_trt.mean() - se_con.mean()) / \
                    (se_trt.std() + se_con.std())
            return res

        return df.apply(
            lambda x: Signal2Noise_method(x[idx_trt], x[idx_con]),
            axis=1,
        )


Rank_Score = {
    'geneList': Rnk_Score.genelist,
    'Signal2Noise': Rnk_Score.Signal2Noise,
}

def mk_gset(idata, orgdb, genelist_sep=';', is_convert_geneid=True):
    """
    Make DataFrame of 'TERM2GENE', 'TERM2NAME' for GSEA

    Required:
    idata:<Path> or <DataFrame>
        if <Path>, input path of a 'csv' file.
        if <DataFrame>, input a DataFrame.
        File must has three columns. One for gene set ids, one for gene set description, one for gene list.
    orgdb:<str> genome database that used in 'clusterProfiler', all available 
        database name can be listed by '.genomeDB.lst_GOfuncR_DB' function, 
        default 'org.Hs.eg.db'.

    Optional:
    genelist_sep:<str> the separators between genes of gene list.    
    """

    if isinstance(idata, str):
        idata = pd.read_csv(idata)
    assert isinstance(idata, pd.DataFrame)

    idata.columns = ['ID', 'Description', 'geneID']
    # make Term id to Gene id DataFrame
    res_TERM2GENE = []
    for row in idata[['ID', 'geneID']].values:
        res_TERM2GENE += [
            [row[0], gname] for gname in row[1].split(';')
        ]
    res_TERM2GENE = pd.DataFrame(
        res_TERM2GENE,
        columns=['gset_id', 'SYMBOL']
    )
    res_TERM2GENE = res_TERM2GENE.drop_duplicates(['gset_id', 'SYMBOL'])
    print('<TERM2GENE>\n', res_TERM2GENE.head())
    # convert gene symbols to gene ids
    if is_convert_geneid:
        symbol_ids_lst = convert_geneSymbols(res_TERM2GENE['SYMBOL'], orgdb=orgdb)
        res_TERM2GENE = pd.merge(res_TERM2GENE, symbol_ids_lst, on='SYMBOL')[
            ['gset_id', 'ENTREZID']]

    # make Term id to Term name DataFrame
    res_TERM2NAME = idata[['ID', 'Description']].copy()
    print('<TERM2NAME>\n', res_TERM2NAME.head())

    res = namedtuple('Res', ['TERM2GENE', 'TERM2NAME'])
    return res(res_TERM2GENE.copy(), res_TERM2NAME.copy())


def mk_gset_from_keggEnrich(kegg_res_path, orgdb):
    """Make DataFrame of 'TERM2GENE', 'TERM2NAME' for GSEA from kegg_enrich result.

    Required:
    kegg_res_path: <Path> path to kegg_enrich result.
    orgdb:<str> genome database that used in 'clusterProfiler', all available 
        database name can be listed by '.genomeDB.lst_GOfuncR_DB' function, 
        default 'org.Hs.eg.db'.
    """
    print('<Making DataFrame of "TERM2GENE", "TERM2NAME" for GSEA from:> {}'.format(
        kegg_res_path))

    kegg_res = pd.read_csv(kegg_res_path)

    # make Term id to Gene id DataFrame
    res_TERM2GENE = []
    for row in kegg_res[['ID', 'geneID']].values:
        res_TERM2GENE += [
            [row[0], gname] for gname in row[1].split(';')
        ]
    res_TERM2GENE = pd.DataFrame(
        res_TERM2GENE,
        columns=['kegg_id', 'SYMBOL']
    )
    print('<TERM2GENE>\n', res_TERM2GENE.head())
    # convert gene symbols to gene ids
    symbol_ids_lst = convert_geneSymbols(res_TERM2GENE['SYMBOL'], orgdb=orgdb)
    res_TERM2GENE = pd.merge(res_TERM2GENE, symbol_ids_lst, on='SYMBOL')[
        ['kegg_id', 'ENTREZID']]

    # make Term id to Term name DataFrame
    res_TERM2NAME = kegg_res[['ID', 'Description']].copy()
    print('<TERM2NAME>\n', res_TERM2NAME.head())

    res = namedtuple('Res', ['TERM2GENE', 'TERM2NAME'])
    return res(res_TERM2GENE.copy(), res_TERM2NAME.copy())


def run_gsea(ifile, odir, TERM2GENE, TERM2NAME, has_header=False, idx_trt=None, idx_con=None, rnk_score_method='Signal2Noise', orgdb='org.Hs.eg.db', nPerm=1000, minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, is_output_plot=True, is_convert_symbol=True, is_convert_symbol_back=True):
    """
    Required:
    ifile: <Path> file (in 'csv' format) path of gene list with expression
        value, first column must store gene symbols, second column should be 
        numerical values for ranking (eg. gene expression values).
    TERM2GENE: <DataFrame> user input annotation of TERM TO GENE mapping, a 
        data.frame of 2 column with term and gene.
    TERM2NAME: <DataFrame> user input of TERM TO NAME mapping, a data.frame of 
        2 column with term and name.
    idx_trt: <str> a list of column index separated by ',' of data in treatment 
        group. Column index example (index starts with 0): '0,1,2,3'.  
    idx_con: <str> a list of column index separated by ',' of data in control 
        group. Column index example (index starts with 0): '0,1,2,3'.  

    Optional:
    orgdb: <str> genome database that used in 'clusterProfiler', all available 
        database name can be listed by '.genomeDB.lst_GOfuncR_DB' function, 
        default 'org.Hs.eg.db'.
    nPerm: <int> number of permutations. Default 1000.
    minGSSize: <int> minimal size of each geneSet for analyzing. Default 10.
    maxGSSize: <int> maximal size of genes annotated for testing. Default 500.
    pvalueCutoff: <float> pvalue cutoff. Default 0.05.
    is_convert_symbol_back: <bool> if convert gene symbols of GO term to 
        original symbols. Default 'True'.
    """
    symbol_ids_lst=None

    # make a new directory if odir not exist
    if not os.path.exists(odir):
        os.mkdir(odir)
        print('<Make a new directory:> "{}" .'.format(odir))
    # check which R is used
    print('<Now using R from:> {}'.format(ro.r('R.home()')))

    path_gsea_obj = os.path.join(odir, 'gsea_res.Rdata')
    path_gene_ann = os.path.join(odir, 'gene_ann.csv')
    path_gsea_res = os.path.join(odir, 'gsea_res.csv')

    print('<Write out "TERM2GENE", "TERM2NAME">')
    path_term2gene = os.path.join(odir, 'term2gene.csv')
    path_term2name = os.path.join(odir, 'term2name.csv')
    TERM2GENE.to_csv(path_term2gene, header=True, index=False)
    TERM2NAME.to_csv(path_term2name, header=True, index=False)

    header = 0 if has_header else None

    gexp = pd.read_csv(ifile, header=header, index_col=0)
    print('<readin expfile:>\n', gexp.head())

    # Calculate rank score
    idx_trt = [int(i) - 1 for i in idx_trt.split(',')]
    idx_con = [int(i) - 1 for i in idx_con.split(',')]
    glist = Rank_Score[rnk_score_method](gexp, idx_trt, idx_con)

    # convert gene symbols to gene ids
    if is_convert_symbol:
        symbol_ids_lst = convert_geneSymbols(glist.index, orgdb=orgdb)
    # convert to StrVector format
    glist.rename('geneList', inplace=True)
    if is_convert_symbol:
        glist = pd.merge(
            symbol_ids_lst, glist,
            left_on=symbol_ids_lst.columns[0],
            right_index=True)
    #gene_ids = symbol_ids_lst.iloc[:, 1].astype('str')

    # sort gene list
    if is_convert_symbol:
        glist.sort_values(
            glist.columns.tolist(),
            axis=0,
            ascending=False,
            inplace=True,
        )
    else:
        glist.sort_values(
                axis=0,
                ascending=False,
                inplace=True,
            )

    with NamedTemporaryFile('w+t') as f:
        # Write out geneList to temp file
        if is_convert_symbol:
            glist = glist.iloc[:, 1:]
            glist.drop_duplicates(glist.columns[0], inplace=True)
        else:
            glist.drop_duplicates(inplace=True)
        print('<Gene list for running GSEA:>\n', glist.head())
        if is_convert_symbol:
            glist.to_csv(f.name, header=None, index=False)
        else:
            glist.to_csv(f.name, header=None)

        # convert Pandas DataFrame to R DataFrame
        gsea_res = ro.r(
            """            
            library(clusterProfiler)
            
            ## convert csv file to geneList format
            d = read.csv('{ifile}', header=F)
            ## assume 1st column is ID
            ## 2nd column is FC

            ## feature 1: numeric vector
            geneList = d[,2]
            ## feature 2: named vector
            names(geneList) = as.character(d[,1])
            ## feature 3: decreasing orde
            geneList = sort(geneList, decreasing = TRUE)
            print(head(geneList))
            
            #data(geneList, package="DOSE")
            gk.res <- GSEA(
                geneList, 
                exponent=1, 
                nPerm={nPerm}, 
                minGSSize={minGSSize},
                maxGSSize={maxGSSize}, 
                pvalueCutoff={pvalueCutoff}, 
                pAdjustMethod="BH",
                TERM2GENE=read.csv('{path_term2gene}', header=T), 
                TERM2NAME=read.csv('{path_term2name}', header=T), 
                verbose=TRUE, 
                seed=FALSE,
                by="fgsea"
            )

            save(gk.res, file='{save_gsea_res}') 
            print('<Result of gseKEGG named "gk.res" has saved in:> "{save_gsea_res}", which can be added to R work envirenment by using "load()".')

            data.frame(gk.res@result)
            """.format(
                ifile=f.name,
                save_gsea_res=path_gsea_obj,
                nPerm=nPerm,
                minGSSize=minGSSize,
                maxGSSize=maxGSSize,
                pvalueCutoff=pvalueCutoff,
                path_term2gene=path_term2gene,
                path_term2name=path_term2name,
            )
        )

    # convert R DataFrame result to Py DataFrame format
    gsea_res = pd.DataFrame(
        np.array(gsea_res).T,
        columns=list(gsea_res.colnames)
    )

    # convert gene symbols of GO term to original symbols
    if is_convert_symbol_back:
        print('<Convert "ENTREZID" to "Gene Names"> ...')
        ids_symbols_dict = dict(symbol_ids_lst.iloc[:, [1, 0]].values)
        gsea_res['core_enrichment'] = gsea_res['core_enrichment'].map(
            lambda x: ';'.join(
                [ids_symbols_dict[i] for i in x.split('/')]
            )
        )
        # Write out gene symbols to gene ids list
        print('<Writing out to: {}>'.format(path_gene_ann))
        symbol_ids_lst.to_csv(path_gene_ann, index=False)

    # Write out
    print('<Writing out to: {}>'.format(path_gsea_res))
    gsea_res.to_csv(path_gsea_res, index=False)

    res = namedtuple('Res', ['gsea_res', 'convert_res'])

    if gsea_res.empty:
        print('<Note:> There is no gene set significant enriched!')
        return res(gsea_res, symbol_ids_lst)

    # plot GSEA result
    if is_output_plot:
        odir_plots = os.path.join(odir, 'gsea_plots')
        if not os.path.exists(odir_plots):
            os.mkdir(odir_plots)
            print('<Make directory for GSEA plots:> "{}" .'.format(odir_plots))
        for i, (gset_id, gset_name) in enumerate(gsea_res[['ID', 'Description']].values):
            path_plots = os.path.join(odir_plots, 'gsea_{}.png'.format(gset_id))
            ro.r(
                """
                library(ggplot2)
                library(clusterProfiler)

                # load gseKEGG result
                load('{load_gsekegg_res}')

                # plot GSEA result
                enrichplot::gseaplot2(
                    gk.res, geneSetID={geneSetID},
                    title='{gset_name}',
                )
                ggsave('{path_plots}')
                """.format(
                    load_gsekegg_res=path_gsea_obj,
                    geneSetID=i + 1,
                    path_plots=path_plots,
                    gset_name=gset_name,
                )
            )
        print('<Total output plots numbers:>', i + 1)

    return res(gsea_res, symbol_ids_lst)

def plot_gsea(gset_info, path_gsea_dir, path_plot):
    """
    gset_info: <iter> first element is gset identity, second element is gset identity type ('ID' or 'Des')
    """
    path_gsea_res = os.path.join(path_gsea_dir, 'gsea_res.csv')
    path_gsea_obj = os.path.join(path_gsea_dir, 'gsea_res.Rdata')

    gset_type2gsea_col_dict = {
        'ID': 'ID',
        'Des': 'Description',
    }

    gset_id, gset_type = gset_info

    gsea_res = pd.read_csv(path_gsea_res)
    gset_id_list = gsea_res[gset_type2gsea_col_dict[gset_type]]
    #print(gset_id_list)

    assert gset_id in list(gset_id_list)
    geneSetID = int(gsea_res[gset_id_list == gset_id].index[0]) + 1
    gset_name = str(gsea_res[gset_id_list == gset_id]['Description'].values[0])

    ro.r(
        """
        library(ggplot2)
        library(clusterProfiler)

        # load gseKEGG result
        load('{load_gsekegg_res}')

        # plot GSEA result
        enrichplot::gseaplot2(
            gk.res, geneSetID={geneSetID},
            title='{gset_name}',
        )
        ggsave('{path_plots}')
        """.format(
            load_gsekegg_res=path_gsea_obj,
            geneSetID=geneSetID,
            path_plots=path_plot,
            gset_name=gset_name,
        )
    )



def gsea_kegg(ifile, odir, idx_trt, idx_con, has_header=False, rnk_score_method='Signal2Noise', orgdb='org.Hs.eg.db', organism='hsa', nPerm=1000, minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, is_convert_symbol_back=True):
    """Run GSEA using KEGG.
    
    Required:
    ifile: <Path> file (in 'csv' format) path of gene list with expression
        value, first column must store gene symbols, second column should be 
        numerical values for ranking (eg. gene expression values).
    odir: <Path> file path to write out all KEGG enrichment result. If not 
        provide, None will be write out. Default 'None'.
    idx_trt: <str> a list of column index separated by ',' of data in treatment 
        group. Column index example (index starts with 0): '0,1,2,3'.  
    idx_con: <str> a list of column index separated by ',' of data in control 
        group. Column index example (index starts with 0): '0,1,2,3'.    
    
    Optional:
    has_header: <bool> if have header or not, default 'False'.
    rnk_score_method: <str> rank score method of object 'Rnk_Score', can be one 
        of 'Signal2Noise, geneList'. Default 'geneList.
    orgdb: <str> genome database that used in 'clusterProfiler', all available 
        database name can be listed by '.genomeDB.lst_GOfuncR_DB' function, 
        default 'org.Hs.eg.db'.
    organism: <str> organism of KEGG. All available organisms can be find in 
        "https://www.genome.jp/kegg/catalog/org_list.html".   
    nPerm: <int> permutation numbers, default 1000.
    minGSSize: <int> minimal size of each geneSet for analyzing, default 10.
    maxGSSize: <int> maximal size of genes annotated for testing, default 500.
    pvalueCutoff: <float> pvalue Cutoff, default 0.05.
    is_convert_symbol_back: <bool> if convert gene symbols of GO term to 
        original symbols. Default 'True'.
    """
    # make a new directory if odir not exist
    if not os.path.exists(odir):
        os.mkdir(odir)
        print('<Make a new directory:> "{}" .'.format(odir))
    # check which R is used
    print('<Now using R from:> {}'.format(ro.r('R.home()')))

    path_gsekegg_obj = os.path.join(odir, 'gsekegg_res.Rdata')
    path_gene_ann = os.path.join(odir, 'gene_ann.csv')
    path_gsekegg_res = os.path.join(odir, 'gsekegg_res.csv')

    header = 0 if has_header else None

    gexp = pd.read_csv(ifile, header=header, index_col=0)
    print('<readin expfile:>\n', gexp.head())

    # Calculate rank score
    idx_trt = [int(i) - 1 for i in idx_trt.split(',')]
    idx_con = [int(i) - 1 for i in idx_con.split(',')]
    glist = Rank_Score[rnk_score_method](gexp, idx_trt, idx_con)

    # convert gene symbols to gene ids
    symbol_ids_lst = convert_geneSymbols(glist.index, orgdb=orgdb)
    # convert to StrVector format
    glist.rename('geneList', inplace=True)
    glist = pd.merge(
        symbol_ids_lst, glist,
        left_on=symbol_ids_lst.columns[0],
        right_index=True)
    #gene_ids = symbol_ids_lst.iloc[:, 1].astype('str')

    # sort gene list
    glist.sort_values(
        glist.columns.tolist(),
        axis=0,
        ascending=False,
        inplace=True,
    )

    with NamedTemporaryFile('w+t') as f:
        # Write out geneList to temp file
        glist = glist.iloc[:, 1:]
        glist.drop_duplicates(glist.columns[0], inplace=True)
        print('<Gene list for running GSEA:>\n', glist.head())
        glist.to_csv(f.name, header=None, index=False)

        gk_res = ro.r(
            """            
            library(clusterProfiler)
            
            ## convert csv file to geneList format
            d = read.csv('{ifile}', header=F)
            ## assume 1st column is ID
            ## 2nd column is FC

            ## feature 1: numeric vector
            geneList = d[,2]
            ## feature 2: named vector
            names(geneList) = as.character(d[,1])
            ## feature 3: decreasing orde
            geneList = sort(geneList, decreasing = TRUE)
            print(head(geneList))
            
            #data(geneList, package="DOSE")
            gk.res <- gseKEGG(
                geneList,
                organism = '{organism}', 
                keyType = "kegg", 
                exponent = 1,
                nPerm = {nPerm}, 
                minGSSize = {minGSSize}, 
                maxGSSize = {maxGSSize},
                pvalueCutoff = {pvalueCutoff}, 
                pAdjustMethod = "BH", 
                verbose = TRUE,
                use_internal_data = FALSE, 
                seed = FALSE, 
                by = "fgsea"
            )

            save(gk.res, file='{save_gsekegg_res}') 
            print('<Result of gseKEGG named "gk.res" has saved in:> "{save_gsekegg_res}", which can be added to R work envirenment by using "load()".')

            data.frame(gk.res@result)
            """.format(
                ifile=f.name,
                save_gsekegg_res=path_gsekegg_obj,
                organism=organism,
                nPerm=nPerm,
                minGSSize=minGSSize,
                maxGSSize=maxGSSize,
                pvalueCutoff=pvalueCutoff,
            )
        )

    # convert R DataFrame result to Py DataFrame format
    gk_res = pd.DataFrame(
        np.array(gk_res).T,
        columns=list(gk_res.colnames)
    )

    # convert gene symbols of GO term to original symbols
    if is_convert_symbol_back:
        print('<Convert "ENTREZID" to "Gene Names"> ...')
        ids_symbols_dict = dict(symbol_ids_lst.iloc[:, [1, 0]].values)
        gk_res['core_enrichment'] = gk_res['core_enrichment'].map(
            lambda x: ';'.join(
                [ids_symbols_dict[i] for i in x.split('/')]
            )
        )

    # Write out gene symbols to gene ids list
    print('<Writing out to: {}>'.format(path_gene_ann))
    symbol_ids_lst.to_csv(path_gene_ann, index=False)
    # Write out
    print('<Writing out to: {}>'.format(path_gsekegg_res))
    gk_res.to_csv(path_gsekegg_res, index=False)

    res = namedtuple('Res', ['gk_res', 'convert_res'])

    if gk_res.empty:
        print('<Note:> There is no gene set significant enriched!')
        return res(gk_res, symbol_ids_lst)

    # plot GSEA result
    odir_plots = os.path.join(odir, 'gsea_plots')
    if not os.path.exists(odir_plots):
        os.mkdir(odir_plots)
        print('<Make directory for GSEA plots:> "{}" .'.format(odir_plots))
    for i, (kegg_id, kegg_name) in enumerate(gk_res[['ID', 'Description']].values):
        path_plots = os.path.join(odir_plots, 'gsea_{}.png'.format(kegg_id))
        ro.r(
            """
            library(ggplot2)
            library(clusterProfiler)

            # load gseKEGG result
            load('{load_gsekegg_res}')

            # plot GSEA result
            enrichplot::gseaplot2(
                gk.res, geneSetID={geneSetID},
                title='{kegg_name}',
            )
            ggsave('{path_plots}')
            """.format(
                load_gsekegg_res=path_gsekegg_obj,
                geneSetID=i + 1,
                path_plots=path_plots,
                kegg_name=kegg_name,
            )
        )
    print('<Total output plots numbers:>', i + 1)

    return res(gk_res, symbol_ids_lst)

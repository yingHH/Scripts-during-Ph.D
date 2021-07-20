# -*- coding: utf-8 -*-
"""
Author: Ying Huang

Date: 2020-03-16 17:26:36
Last Modified by: Ying Huang
Last Modified time: 2020-03-16 17:26:36

Description: 
Find DEGs (differential expression genes) by using R package DEseq2
"""
import os
import pandas as pd
import rpy2.robjects as ro


def count_degs(countfile, odir, trt_cols, ctl_cols):
    """
    count DEGs between treat and control.

    Required:
    countfile: <Path> reads count matrix file (csv) path. First column must be gene ids,
        the other columns are reads count values. First row must be column names.
    odir: <Path> result output directory.
    trt_cols: <tuple> column names for treatment samples.
    ctl_cols: <tuple> column names for control samples (starts with 0).
    """
    if not os.path.exists(odir):
        os.mkdir(odir)
    # make coldata file based on countfile
    colfile = os.path.join(odir, 'coldata.csv')

    condition = pd.DataFrame(
        [[i,'treated'] for i in trt_cols] + [[i,'untreated'] for i in ctl_cols],
        columns=['samples', 'condition']
    )
    print('<Write coldata to> :', colfile)
    condition.to_csv(colfile, index=False)

    # run 'DESeq2'
    ofile = os.path.join(odir, 'DEGs.deseq2.csv')
    save_degs = os.path.join(odir, 'DEGs.deseq2.Rdata')

    print('<Run DESeq2> ...')

    # check R settings
    ro.r(
        """
        r_home = R.home()
        r_libPath = .libPaths()
        """
    )
    print(
        """
        > using R_HOME: {r_home}
        > using R_libPath: {r_libPath}
        """.format(
            r_home=ro.globalenv['r_home'].r_repr(),
            r_libPath=ro.globalenv['r_libPath'].r_repr(),
        )
    )

    ro.r(
        """
        # run DESeq2 with 'BGC vs SGC'
        # ref "https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual"
        # ref "https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input"
        library(DESeq2)

        # read 'countData' and 'coldata'
        countData <- as.matrix(read.csv('{countfile}', row.names=1))
        print(head(countData), quote = TRUE, row.names = FALSE)
        countData_head = head(countData)

        colData <- read.csv('{colfile}', row.names=1)
        print(head(colData), quote = TRUE, row.names = FALSE)
        colData_head = head(colData)

        # out name
        oname <- '{ofile}'

        # Check all sample IDs in colData are also in CountData and match their orders
        print(all(rownames(colData) %in% colnames(countData)), quote = TRUE, row.names = FALSE)

        # Create a DESeqDataSet
        dds <- DESeqDataSetFromMatrix(countData = countData,
                                    colData = colData,
                                    design = ~ condition)
        print(dds, row.names = FALSE)

        # Run DESeq2
        dds <- DESeq(dds)
        res <- results(dds)
        resOrdered <- res[order(res$padj), ]

        # write out
        cat("<Write result to> :", oname)
        write.csv(resOrdered, oname)

        save(dds, file='{save_degs}')
        """.format(
            countfile = countfile,
            colfile=colfile,
            ofile=ofile,
            save_degs=save_degs,
        )
    )

    print(
       """
       <Write result to> : {ofile}
       <Save DESeq2 'dds' to> : {save_degs}
       """.format(
           ofile=ofile,
           save_degs=save_degs,
       )
    )
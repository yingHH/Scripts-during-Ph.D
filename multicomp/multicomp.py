# -*- coding: utf-8 -*-
"""
Created on 2021-05-19 14:32:31
Last Modified on 2021-05-19 14:32:31

Running multiple comparisons by using R package 'agricolae' 

@Author: Ying Huang
"""
import os
import pandas as pd
import rpy2.robjects as ro
from collections import namedtuple


def mcmp(ifile, value, grp, oname, alpha=0.05, padj='bonferroni'):

    Res = namedtuple('Res', ['pval', 'mcomp'])

    r_script = """
        library(agricolae)

        data <- read.csv('{ifile}')
        print(head(data))

        model <- aov(
            formula = {value} ~ {grp},
            data = data
        )
        print(summary(model))
        # writ anova result to txt
        sink("{oname}.anova_res.txt")
        print(summary(model))
        sink()

        # get and write multi-compare pvalue
        cmp.syf <- LSD.test(model, '{grp}', alpha = {alpha}, p.adj = "{padj}", group = F)
        write.csv(cmp.syf$comparison, file = '{oname}.pval_res.csv')
        # get and write multi-compare result
        res.syf <- LSD.test(model, '{grp}', alpha = {alpha}, p.adj = "{padj}")
        write.csv(res.syf$groups, file = '{oname}.multicomp_res.csv')
        """.format(
        ifile=ifile,
        value=value,
        grp=grp,
        alpha=alpha,
        padj=padj,
        oname=oname,
    )
    #print(r_script)

    ro.r(
        r_script
    )

    print(
        '> Results are: \n [1] "{oname}.pval_res.csv"\n [2] "{oname}.multicomp_res.csv"\n [3] "{oname}.anova_res.txt".'.format(
            oname=oname)
    )

    pval = pd.read_csv('{oname}.pval_res.csv'.format(oname=oname), index_col=0)
    mcomp = pd.read_csv(
        '{oname}.multicomp_res.csv'.format(oname=oname), index_col=0)

    return Res(pval, mcomp)

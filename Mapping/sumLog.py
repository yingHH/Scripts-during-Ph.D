# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 16:07:20 2018

@author: yinghuang
"""
import fire
import os
import pandas as pd

from getLogInfo import get_log_info


def main(idir, ofilename='summary.log.csv'):
    sum_logFiles(idir, ofilename)

def sum_logFiles(idir, ofilename):
    logDf = sum_log(idir)
    logDf.to_csv(os.path.join(idir, ofilename))
    print("> Statistic all log files and write to :", ofilename)


def sum_log(idir):
    fileLst = lst_log_files(idir)
    sumLog = statDfLog(fileLst)
    sumLog.count_dup_rate()
    return sumLog.sumDf
    #print(sumLog.sumDf)


def lst_log_files(idir):
    lst = [os.path.join(os.path.abspath(idir), f)
           for f in os.listdir(idir)
           if f.endswith('.log')]
    return tuple(lst)


class ConcatDfLog:
    def __init__(self, fileLst):
        self.flst = fileLst
        self.sumDf = self._sum()

    def _get_sampName(self, ifileName):
        return os.path.basename(ifileName).rstrip('.log')

    def _normalize(self, dfLog, sampName):
        #print('df [_normalize] :\n', dfLog.head())
        df = dfLog.set_index([0])
        df.rename(columns = {1: sampName}, inplace=True)
        return df

    def _concat_df(self, df1, df2):
        return pd.concat([df1, df2], axis=1)

    def _sum(self):
        sumDf = pd.DataFrame()
        for f in self.flst:
            try:
                print("> [log file] :", f)
                sampName = self._get_sampName(f)
                dfLog = get_log_info(f)
                #print("[_sum][get_log_info] :", dfLog.head())
                #print('get_log_info [out]:', dfLog)
                nDfLog = self._normalize(dfLog, sampName)
                sumDf = self._concat_df(sumDf, nDfLog)
            except:
                print("> [err][failed to statistic log file] :", f)
                continue
        return sumDf.T


class statDfLog(ConcatDfLog):
    def count_dup_rate(self):
        #print("[statDfLog] :\n", self.sumDf)
        self.sumDf['DUPLICATE Reads'] = \
        self.sumDf['READ Reads'] - self.sumDf['WRITTEN Reads']

        self.sumDf['DUPLICATE Rate'] = \
        self.sumDf['DUPLICATE Reads'] / self.sumDf['READ Reads']
        self.sumDf['DUPLICATE Rate'] = self.sumDf['DUPLICATE Rate'].map(
            lambda x: str(format(x * 100, '0.2f')) + '%')


if __name__ == '__main__':
    fire.Fire(main)

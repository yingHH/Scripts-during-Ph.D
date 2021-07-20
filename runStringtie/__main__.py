#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 17:13:40 2018

@author: yinghuang
"""
from concurrent import futures
import fire
import multiprocessing
import os
import pandas as pd


PWD = os.path.split(os.path.realpath(__file__))[0]


def main(idir, odir, gff3, oname='reads_count', nproc=1):
    if not os.path.isdir(odir):
        os.mkdir(odir)

    run_stringties(idir, odir, gff3, nproc)
    ilst = os.path.join(odir, 'samp_list.txt')
    print("> Running 'prepDE.py' ...")
    run_prepDE(ilst, odir, oname)
    print("> Done.")


def run_stringties(idir, odir, gff3, nproc):
    sampLst = []

    safeProcNum = int(multiprocessing.cpu_count() * 0.9)
    nproc = min(safeProcNum, nproc)
    isCancel = False
    print("> Running 'stringtie's ...")
    with futures.ThreadPoolExecutor(max_workers=nproc) as executor:
        toDo = []
        for ibam in lst_bam(idir):
            # get output file name base on name of bam file
            sampName = os.path.basename(ibam).split('.')[0]
            ogtf = os.path.join(odir, sampName + '.abund.gtf')
            otxt = os.path.join(odir, sampName + '.abund.txt')
            sampLst.append((sampName, ogtf))

            future = executor.submit(run_stringtie, ibam, gff3, ogtf, otxt)
            toDo.append(future)
            msg = 'Scheduled for {}: {}'
            print(msg.format(ibam, future))

        try:
            results = []
            for future in futures.as_completed(toDo):
                res = future.result()
                msg = '{} result: {!r}'
                print(msg.format(future, res))
                results.append(res)
        except KeyboardInterrupt:
            isCancel = True
            for future in futures:
                future.cancel()

        if isCancel:
            executor.shutdown()

    # write 'sampLst file' for 'prepDE.py' input
    print("> Writing 'samp_list.txt' file ...")
    opathSampLst = os.path.join(odir, 'samp_list.txt')
    pd.DataFrame(sampLst).to_csv(opathSampLst, header=None, index=False, sep="\t")

    return len(results)


def sampName(ipath):
    return os.path.basename(ipath).split('.')[0]

def lst_bam(idir):
    for f in os.listdir(idir):
        if f[-4:] == '.bam':
            yield os.path.join(idir, f)


def run_stringtie(ifile, gff3, ogtf, otxt):
    cmd = 'bash {} {} {} {} {}'
    script = os.path.join(PWD, 'stringtie.pip.sh')
    signal = os.system(cmd.format(script, ifile, gff3, ogtf, otxt))
    return signal


def run_prepDE(ilst, odir, oname):
    cmd = 'bash {} {} {} {}'
    script = os.path.join(PWD, 'prepDE.pip.sh')
    signal = os.system(cmd.format(script, ilst, odir, oname))
    return signal


if __name__ == '__main__':
    fire.Fire(main)

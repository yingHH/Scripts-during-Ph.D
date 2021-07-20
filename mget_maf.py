# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 22:57:31 2021

@author: Ying
"""
from Bio import AlignIO
import requests
import gzip
import os
import pandas as pd
from concurrent import futures
import multiprocessing


def wget_maf(filename, odir, url = "http://hgdownload.soe.ucsc.edu/goldenPath/galGal6/multiz77way/maf/"):
    print('> Download file: ' + filename)
    rfile = requests.get(url + filename, allow_redirects=True)
    
    print('> Ungzip file: ' + filename)
    ug_file = gzip.decompress(rfile.content)
    
    ofile = os.path.join(odir, filename.rstrip('\.gz'))
    with open(ofile, 'wb') as f:
        f.write(ug_file)
        
    return ofile
        
def make_mafidx(maf, seq_name):
    idx_name = maf + 'index'
    print('> Make mafindex with parameters: "{}", "{}", "{}"'.format(idx_name, maf, seq_name))
    AlignIO.MafIO.MafIndex(idx_name, maf, seq_name)
    
def get_maf(filename, odir, sp, url = "http://hgdownload.soe.ucsc.edu/goldenPath/galGal6/multiz77way/maf/"):
    ofile = wget_maf(filename, odir, url)
    maf = ofile
    seq_name = '{}.{}'.format(sp, filename[:-7])
    make_mafidx(maf, seq_name)
    
def mget(ilist, odir, max_proc=5):
    todo_df = pd.read_csv(ilist)
    
    isCancel = False
    with futures.ThreadPoolExecutor(max_workers=max_proc) as executor:
        toDo = []
        for filename, sp, url in todo_df.values.tolist():
            future = executor.submit(get_maf, filename, odir, sp, url)
            toDo.append(future)
            msg = 'Scheduled for {}: {}'
            #print(msg.format(filename, future))

        try:
            results = []
            for future in futures.as_completed(toDo):
                res = future.result()
                msg = '{} result: {!r}'
                #print(msg.format(future, res))
                results.append(res)
        except KeyboardInterrupt:
            isCancel = True
            for future in futures:
                future.cancel()

        if isCancel:
            executor.shutdown()
    return len(results)
    
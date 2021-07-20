# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 21:27:21 2018

@author: yinghuang
"""
from concurrent import futures
import multiprocessing
import os
import pandas as pd

from pipABC import RunPip


class CtrlPip:

    def __init__(self, pipName, nproc, infoCsv, refgenome, idir, odir):
        """
        pipName:str
        nproc:int
        """
        self.pipName = pipName
        safeProcNum = int(multiprocessing.cpu_count() * 0.9)
        self.nproc = min(nproc, safeProcNum)
        self.fileNames = pd.read_csv(infoCsv).iloc[:, [0, 1]].values.tolist()
        self.refgenome = refgenome
        self.idir = idir
        self.odir = odir
        ifiles = [(os.path.join(self.idir, fileName[0]), os.path.join(self.idir, fileName[1]))
                       for fileName in self.fileNames]
        self.ifiles = [(os.path.abspath(f[0]), os.path.abspath(f[1])) for f in ifiles]
        self.all_seq_type = set()

    def do_one_pip(self, pipName, ifiles, odir, refgenome):
        self.pipRes = RunPip(pipName, ifiles, odir, refgenome)
        self.all_seq_type.update(self.pipRes.seq_type)
        return str(self.pipRes)

    def do_pips(self):
        isCancel = False
        with futures.ThreadPoolExecutor(max_workers=self.nproc) as executor:
            toDo = []
            for fs in self.ifiles:
                future = executor.submit(self.do_one_pip, self.pipName, fs, self.odir, self.refgenome)
                toDo.append(future)
                msg = 'Scheduled for {}: {}'
                print(msg.format(fs, future))

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
        return len(results)

    def __repr__(self):
        res = self.do_pips()
        return "> Total exe number {}".format(res)

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 21:27:21 2018

@author: yinghuang
"""
from concurrent import futures
import multiprocessing
import os

from pipABC import RunPip


class CtrlPip:

    def __init__(self, shName, nproc, idir, odir):
        """
        shName:str
        nproc:int
        """
        self.shName = shName
        safeProcNum = int(multiprocessing.cpu_count() * 0.9)
        self.nproc = min(nproc, safeProcNum)
        self.idir = idir
        self.odir = odir
        self.ifiles = [os.path.join(self.idir, fileName)
                       for fileName in os.listdir(self.idir)
                       if fileName[-3:] == 'bam']

    def do_one_pip(self, shName, ifiles, odir):
        self.pipRes = RunPip(shName, ifiles, odir)
        return str(self.pipRes)

    def do_pips(self):
        isCancel = False
        with futures.ThreadPoolExecutor(max_workers=self.nproc) as executor:
            toDo = []
            for file in self.ifiles:
                future = executor.submit(self.do_one_pip, self.shName, file, self.odir)
                toDo.append(future)
                msg = 'Scheduled for {}: {}'
                print(msg.format(file, future))

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
        return "Total exe number {}".format(res)

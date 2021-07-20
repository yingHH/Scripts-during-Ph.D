#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 21:27:21 2018

@author: yinghuang
"""
import os
import fire

from runpips import CtrlPip


def main(shName, nproc, idir, odir):
    print('starting ...')
    if not os.path.exists(idir):
        raise Exception("Path '{}' not exist.".format(idir))
    elif not os.path.exists(odir):
        print("Make directory '{}'.".format(odir))
        os.mkdir(odir)

    res = CtrlPip(shName, nproc, idir, odir)
    print(str(res))


if __name__ == '__main__':
    fire.Fire(main)

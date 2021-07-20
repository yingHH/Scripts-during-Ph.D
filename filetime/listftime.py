# -*- coding: utf-8 -*-
"""
Created on 2019-09-08 19:44:06
Last Modified on 2019-09-08 19:44:06

List create time, access time and modified time of a set of files in a directory.

@Author: Ying Huang
"""
import os
import time
import pandas as pd


def list_file_time(idir, time_type='ctime', ofile=None):
    """Get create, access, modified time of a directory.
    
    Required: 
    idir: <Path> path of a file or directory.
    ofile: <Path> path to write out result.

    Optional:
    time_type: <str> one of 'ctime', 'atime', 'mtime'. In detail, 'ctime'           (create time), 'atime'(access time), 'mtime' (modified time), default       'ctime'."""

    assert time_type in ['ctime', 'atime', 'mtime']
    assert os.path.isdir(idir)

    time_methods = {
        'ctime': os.path.getctime,
        'atime': os.path.getatime,
        'mtime': os.path.getmtime,
    }

    list_of_dfiles = os.listdir(idir)
    ftimes_of_dfiles = []

    for ifile in list_of_dfiles:
        ifile = os.path.join(idir, ifile)
        ftime = time.localtime(
            time_methods[time_type](ifile)
        )
        fmt_ftime = time.strftime('%Y-%m-%d', ftime)

        ftimes_of_dfiles.append(fmt_ftime)

    res = pd.DataFrame([ftimes_of_dfiles, list_of_dfiles]).T
    res.columns = [time_type, 'files']

    return res.sort_values(time_type).copy()

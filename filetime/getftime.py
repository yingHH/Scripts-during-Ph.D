# -*- coding: utf-8 -*-
"""
Created on 2019-09-08 16:54:21
Last Modified on 2019-09-08 16:54:21

Get create time, access time and modified time of a file/directory.

@Author: Ying Huang
"""
import os
import time


def get_file_time(ifile, time_type='ctime'):
    """Get create, access, modified time of a file/directory.
    
    Required: 
    ifile: <Path> path of a file or directory.

    Optional:
    time_type: <str> one of 'ctime', 'atime', 'mtime'. In detail, 'ctime'           (create time), 'atime'(access time), 'mtime' (modified time), default       'ctime'."""

    assert time_type in ['ctime', 'atime', 'mtime']
    assert os.path.exists(ifile)

    time_methods = {
        'ctime': os.path.getctime,
        'atime': os.path.getatime,
        'mtime': os.path.getmtime,
    }

    ftime = time.localtime(
        time_methods[time_type](ifile)
    )

    return time.strftime('%Y-%m-%d', ftime)



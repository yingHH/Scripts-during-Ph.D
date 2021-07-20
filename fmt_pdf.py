# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 09:48:19 2021

@author: Administrator
"""

import os
import re


def fmt_pdf(idir):
    for (root , dirs, files) in os.walk(idir):
        for f in files:
            absf = os.path.join(root, f)
            if absf.endswith('.pdf'):
                newf = re.sub('(\.pdf)+', '.pdf', r'{}'.format(absf))
                if (absf != newf) and (not os.path.exists(newf)):
                    print(">Rename: '{}' -> '{}'".format(absf, newf))
                    os.rename(absf, newf)
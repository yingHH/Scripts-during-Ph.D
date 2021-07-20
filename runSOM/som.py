# -*- coding: utf-8 -*-
"""
Created on 2019-10-22 09:20:05
Last Modified on 2019-10-22 09:20:05

A core

@Author: Ying Huang
"""
from collections import namedtuple
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from tempfile import TemporaryDirectory
import os


def kohonen_r():


    with TemporaryDirectory() as dirname:
        
        som_res = ro.r(
            """
            library('kohonen')

            idata = read.csv('{ifile}', sep = ",", header=T, row.names=1)

            idata_z = scale(idata)
            print(head(idata_z))

            idata_som = supersom(
                nba_z, 
                grid=somgrid(
                    {grid_x}, {grid_y}, 
                    "hexagonal", toroidal=T,
                ), 
                core={nproc}
            )

            
            
            """
        )

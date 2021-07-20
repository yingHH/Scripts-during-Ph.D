# -*- coding: utf-8 -*-
"""
Created on Wed May 23 16:01:47 2018

@author: yinghuang
"""
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

color_lst = ("#BD0026", "#E31A1C", "#E31A1C", "#E31A1C", "#FC4E2A",
             "#1A9850", "#006837", "#A6D96A", "#A6D96A", "#FEB24C", 
             "#FEB24C", "#FFFF33", "#35978F", "#807DBA", "#B2182B", 
             "#737373", "#969696", "#FFFFFF", "#FFFFFF", "#4169e1")

# https://matplotlib.org/2.0.1/users/legend_guide.html?highlight=plt%20legend
for i, color in enumerate(color_lst):
    plt.plot(0,
            i, 
            'o',
            c=color,
            markersize=10,
            label=color)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
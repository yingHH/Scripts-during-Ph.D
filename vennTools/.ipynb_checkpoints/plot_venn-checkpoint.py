# -*- coding: utf-8 -*-
"""
Author: Ying Huang
Date: 2019-11-10 17:06:03
Last Modified by: Ying Huang
Last Modified time: 2019-11-10 17:06:03
Descriptipn: plot venn2, venn3
"""
# self packges
from .count_venn import count_venn2, count_venn3

def plot_venn2(isets, labels=('grpA', 'grpB'), colors= ('#e74c3c', '#00755E'), linestyle='-', linewidth=1, linecolor="grey", ofile=None, figsize=(6,6), dpi=300):
    """
    Plot 2 venn

    Required:
    isets:<tuple> a tuple contains 2 sets for doing venn

    Optional:
    labels:<tuple> labels of 2 sets. Default:('grpA', 'grpB')
    colors:<tuple> colors of 2 sets. Default:('#e74c3c', '#00755E')
    linestyle:<str> boundry style of 2 sets. Default:'-'
    linewidth:<float> bound of boundry. Default:1
    linecolor:<str> color of boundry. Default:'grey'
    ofile:<str> output file name. Default:None
    figsize:<tuple> figure size. Defalut:(6,6)
    dpi:<int> Default:300
    """

    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2, venn2_circles


    plt.figure(figsize=figsize)

    subsets = count_venn2(isets)

    v=venn2(subsets = subsets, set_labels = labels)
    # set color in each patch
    for label, color in zip(('10', '01'), colors):
        v.get_patch_by_id(label).set_color(color)

    c = venn2_circles(
        subsets = subsets, linestyle=linestyle, linewidth=linewidth, color=linecolor
    )

    if ofile is not None:
        plt.savefig(ofile, dpi=dpi)

    res = plt.gca().figure
    plt.close()

    return res



def plot_venn3(isets, labels=('grpA', 'grpB', 'grpC'), colors=('#e74c3c', '#3498db', '#00755E'), linestyle='-', linewidth=1, linecolor="grey", ofile=None, figsize=(6,6), dpi=300):
    """
    Plot 3 venn

    Required:
    isets:<tuple> a tuple contains 3 sets for doing venn

    Optional:
    labels:<tuple> labels of 3 sets. Default:('grpA', 'grpB', 'grpC')
    colors:<tuple> colors of 3 sets. Default:('#e74c3c', '#3498db', '#00755E')
    linestyle:<str> boundry style of 3 sets. Default:'-'
    linewidth:<float> bound of boundry. Default:1
    linecolor:<str> color of boundry. Default:'grey'
    ofile:<str> output file name. Default:None
    figsize:<tuple> figure size. Defalut:(6,6)
    dpi:<int> Default:300
    """

    import matplotlib.pyplot as plt
    from matplotlib_venn import venn3, venn3_circles
    
    plt.figure(figsize=figsize)

    subsets = count_venn3(isets)

    v=venn3(
        subsets = subsets, set_labels = labels,
    )
    # set color in each patch
    v100 = subsets[0]
    v010 = subsets[1]
    v001 = subsets[3]
    for label, color, count in zip(('100', '010', '001'), colors, [v100, v010, v001]):
        if count > 0:
            v.get_patch_by_id(label).set_color(color)

    c = venn3_circles(
        subsets = subsets, linestyle=linestyle, linewidth=linewidth, color=linecolor
    )

    if ofile is not None:
        plt.savefig(ofile, dpi=dpi)

    res = plt.gca().figure
    plt.close()

    return res
    
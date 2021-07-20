- Here is a example for venn plot.

```py
%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

plt.figure(figsize=[5, 5])

ph = set(crc_heatmap_srt[crc_heatmap_srt[[n for n in crc_heatmap_srt.columns if n.startswith('ph')]].sum(axis=1) > 0].index.tolist())
len(ph)
po = set(crc_heatmap_srt[crc_heatmap_srt[[n for n in crc_heatmap_srt.columns if n.startswith('po')]].sum(axis=1) > 0].index.tolist())
len(po)
mm = set(crc_heatmap_srt[crc_heatmap_srt[[n for n in crc_heatmap_srt.columns if n.startswith('mm')]].sum(axis=1) > 0].index.tolist())
len(mm)

def count_venn3(isets):
    assert len(isets) == 3

    a, b, c = isets

    aINb, aINc, bINc = [i & j for i, j in zip([a, a, b], [b, c, c])]
    abc = a & b & c

    aNum, bNum, cNum, abNum, acNum, bcNum, abcNum = \
    [len(i) for i in [a, b, c, aINb, aINc, bINc, abc]]

    v111 = abcNum
    v110, v101, v011 = [intersect - abcNum for intersect in [abNum, acNum, bcNum]]
    v100, v010, v001 = [whole - interI - interJ - v111
                        for whole, interI, interJ in
                        zip([aNum, bNum, cNum],
                            [v110, v011, v101],
                            [v101, v110, v011])]
   
    return v100, v010, v011, v001, v101, v011, v111

sets = count_venn3([ph, po, mm])
sets

v=venn3(subsets = sets, set_labels = ('ph', 'po', 'mm'))
# set color in each patch
for label, color in zip(('100', '010', '001'), ( '#e74c3c', '#3498db', '#00755E')):
    v.get_patch_by_id(label).set_color(color)

c=venn3_circles(subsets = sets, linestyle='dashed', linewidth=1, color="grey")
plt.savefig('./CRC_all_gephi/mm_all_crcmap.venn.png', dpi=300)
```
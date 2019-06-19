#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

hetero = [428271,359781,344053,355108,361489,350074,346966,688626,
          1360789,2207694,2211361,2294693,1816568,1818054,1787033,1365294,1143449,1396429,1378405,1423698]
homo = [29062,86747,87716,69518,64646,67903,65760, 46967,
        1422983,1235197,1186560,1207352,1194830,1199045,1149691,529898,567360,1446508,560791,56959]

ticks = ["SB6536","SB7462","SB8055","SB6573","SB6815","SB2","SB10", "38741",
         "LIB18989","LIB18990","LIB18991","LIB18992","LIB18993","LIB18994","LIB18995","LIB21974","LIB21975","LIB21977","LIB22032","LIB23764"]


fig = plt.figure(1, figsize=(8,4), dpi=300)

bar_width = 0.33
bin_coord = np.arange(len(ticks))
plt.bar(bin_coord, hetero, width=bar_width)
plt.bar(bin_coord, hetero, width=bar_width, color="blue", label="hetero")
plt.bar(bin_coord + bar_width, homo, width=bar_width, color="green", label="homo")
plt.ylabel('Variants', fontweight='bold')
plt.xlabel('Sample', fontweight='bold')
plt.xticks([coord + bar_width for coord in range(len(bin_coord))], ticks, rotation=45)
plt.legend
plt.legend()
def tick_formater(x, pos):
    return '%.1f M' % (x*1e-6)


formatter = FuncFormatter(tick_formater)
plt.subplot().yaxis.set_major_formatter(formatter)
plt.savefig("variant_zygoty.bff_domfer_epol.png", bbox_inches='tight')

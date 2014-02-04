#!/usr/bin/env python
"""
plot the total energy and mcc value in one plot with two axes
"""


import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

ifn = 'mcc.out'
mcc = pd.read_csv(ifn)


fig, ax1 = plt.subplots()
s1 = mcc['energy']
ax1.plot(s1, 'b')
ax1.set_xlabel('step')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('total energy', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
s2 = mcc['mcc']
ax2.plot(s2, 'r')
ax2.set_ylabel('mcc', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')

fig.savefig('ener_mcc.pdf')

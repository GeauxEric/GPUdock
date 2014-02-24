#!/usr/bin/env python
"""
scatter plot of the total energy and mcc value in one plot with two axes
"""


import pandas as pd
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

ifn = 'baye_mcc.out'
data = pd.read_csv(ifn)
mcc = data['mcc']
ener = data['energy']


fig, ax1 = plt.subplots()
ax1.scatter(mcc, ener, alpha=0.5)
ax1.set_xlabel('mcc')
# ax1.set_ylabel('total energy')
ax1.set_ylabel('probability difference')

fig.savefig('scatter_ener_mcc.pdf')

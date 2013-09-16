"""
lst_hist.py:     losd a list and draw the histogram
"""
import numpy as np
from pandas import Series
import sys

data_file = sys.argv[1]
bins_num = 100 if len(sys.argv) == 2 else int(sys.argv[2])

Series(np.loadtxt(data_file)).hist(bins=bins_num)

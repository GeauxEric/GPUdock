import pandas as pd
from scipy.stats.stats import pearsonr

ifn = 'mcc.out'
data = pd.read_csv(ifn)
mcc = data['mcc']
ener = data['energy']

print pearsonr(mcc, ener)

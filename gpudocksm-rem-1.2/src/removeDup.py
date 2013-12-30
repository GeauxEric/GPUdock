import sys, os
import subprocess
import glob
import pandas
import numpy as np
import matplotlib

import matplotlib.pyplot as plt

df = pandas.DataFrame

def plotTotalEne(dt, pdf_path):
    """line plot of the total energy
    """
    ofn = pdf_path

    ener = dt['total'].values

    fig, ax = plt.subplots()
    ax.plot(ener)
    ax.set_xlabel('step', fontsize=14)
    ax.set_ylabel('total energy', fontsize=14)

    print "line plot of total energy\t", ofn
    plt.savefig(ofn)
    
def getHdfPath(pattern=' '):
    """get the absolute path of the hdf file
    
    Arguments:
    - `pattern`: regex of the hdf file path
    """
    # pattern = './output_*/a_*'
    pattern = pattern
    outputs = glob.glob(pattern)
    hdf_path = outputs[0]  # in testing, only the first hdf file tested
    hdf_path = os.path.abspath(hdf_path)
    
    return hdf_path
    
def main():
    pattern = './output_*/a_*.h5'
    hdf_path = getHdfPath(pattern)
    
    base = os.path.splitext(hdf_path)[0]  # in testing, only the first hdf file tested

    csv_path = base + '.csv'
    total_ener_path = base + '_total.csv'
    pdf_path = base + '.pdf'

    print "pandas loading \t\t\t", csv_path
    dt = pandas.read_csv(csv_path)

    print "removing duplicates ..."
    dt = dt.drop_duplicates(cols=['total', 'vdw', 'ele', 'pmf', 'psp', 'hdb', 'hpc', 'kde', 'lhm', 'dst'])

    print "total enregy\t\t\t", total_ener_path
    total = df(dt['total'])
    total.to_csv(total_ener_path, index=False)

    plotTotalEne(dt, pdf_path)


    
if __name__ == "__main__":
    main()

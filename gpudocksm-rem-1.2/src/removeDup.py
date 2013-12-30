import sys, os
import subprocess
import glob
import pandas
import numpy as np
import matplotlib

import matplotlib.pyplot as plt

df = pandas.DataFrame

def plotTotalEne(dt, hdf_path):
    """line plot of the total energy
    """
    base = os.path.splitext(hdf_path)[0] 
    ofn = base + '.pdf'

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
    

def loadAndRemoveDup(hdf_path):

    base = os.path.splitext(hdf_path)[0] 
    analysis_path = base + '.csv'

    print "pandas loading \t\t\t", analysis_path

    dt = pandas.read_csv(analysis_path)

    print "removing duplicates ..."
    dt = dt.drop_duplicates(cols=['total'])
    
    return dt

def printTotal(dt, hdf_path):
    """print the total energy data to csv 
    """
    base = os.path.splitext(hdf_path)[0] 
    ofn = base + '_total.csv'
    print "total enregy\t\t\t", ofn
    total = df(dt['total'])
    total.to_csv(ofn, index=False)

if __name__ == "__main__":
    pattern = './output_*/a_*.h5'
    hdf_path = getHdfPath(pattern)
    dt = loadAndRemoveDup(hdf_path)
    printTotal(dt, hdf_path)
    plotTotalEne(dt, hdf_path)



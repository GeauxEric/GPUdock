import sys, os
import subprocess
import glob
import pandas
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


df = pandas.DataFrame
pd = pandas

def plotTotalEne(dt, pdf_path):
    """line plot of the total energy
    """
    ofn = pdf_path

    ener = dt['total'].values
    # ener = dt['total'].values

    fig, ax = plt.subplots()
    ax.plot(ener)
    ax.set_xlabel('step', fontsize=14)
    # ax.set_ylabel('total energy', fontsize=14)
    ax.set_ylabel('total energy', fontsize=14)

    print "line plot of total energy\t", ofn
    plt.savefig(ofn)
    
def getCsvPaths(pattern=' '):
    """get the absolute paths of the csv file
    
    Arguments:
    - `pattern`: regex of the hdf file path
    """
    pattern = pattern
    outputs = sorted(glob.glob(pattern))  # sorted by name
    csv_paths = []
    for relative_path in outputs:
        csv_paths.append(os.path.abspath(relative_path))
    
    return csv_paths
    
def loadTotalEner(csv_paths):
    """load all the csv files and extrct the total energy
    """
    total_dt = df()
    for path in csv_paths:
        dt = pd.read_csv(path)
        total_dt = total_dt.append(dt)
        os.remove(path)

    print "removing duplicates ..."
    total_dt = total_dt.drop_duplicates(cols=['total', 'vdw', 'ele', 'pmf', 'psp', 'hdb', 'hpc', 'kde', 'lhm', 'dst'])
    
    total_ener = df(total_dt['total'])
    # total_ener = df(total_dt['total'])

    return total_ener



def main():
    dir_path = sys.argv[1]
    base_name = sys.argv[2]

    pattern = dir_path + '/a_*.csv'
    csv_paths = getCsvPaths(pattern)
    
    total_ener = loadTotalEner(csv_paths)

    pdf_path = dir_path + '/' + base_name + '.pdf'

    plotTotalEne(total_ener, pdf_path)
    
if __name__ == "__main__":
    main()

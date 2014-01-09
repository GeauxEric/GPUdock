import sys, os
import subprocess
import glob

import pandas as pd
import numpy as np

import plot_analysis

df = pd.DataFrame

    
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

    line_path = dir_path + '/' + base_name + '_line.pdf'
    hist_path = dir_path + '/' + base_name + '_hist.pdf'
    values = total_ener['total'].values

    plot_analysis.plotTotalEne(values, line_path)
    plot_analysis.histogram(values, hist_path, normed=0)
    
if __name__ == "__main__":
    main()






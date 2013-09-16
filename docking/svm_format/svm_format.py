import pandas as pd
import numpy as np
import pandas.DataFrame as df
import suffix

def Concat(concat_lst):
    """load a list of file into df, concat them into one df"""
    all_data = []
    for file in concat_lst:
        all_data.append(pd.read_csv(file))
    data = pd.concat(all_data)
    return data

def main():
    """
    read multiple xxx-DecoyRepEner.txt files, concat them
    convert to the format that libsvm accept
    """
    import os
    fullnames = []
    for (thisDirLevel, subsHere, filesHere) in os.walk(rootdir):
        for filename in filesHere:
            if filename.endswith(extd):
                fullnames.append(os.path.abspath(filename))

    paths = suffix.getPaths(rootdir, extd)



    for fullname in fullnames:
        print fullname
    
    


import pandas as pd
from scipy.stats.stats import pearsonr
from readH5 import loadTrack


def pearsonrMccEner(rep_num, h5_path_regx):
    """
    read the hdf5 file and calculate the pearsonr correlation coefficient between mcc and energy
    """
    track = loadTrack(rep_num, h5_path_regx)

    ener, mcc = [], []

    for i in track:
        ener.append(i[1])
        mcc.append(i[2])
        
    return pearsonr(mcc, ener)
    
if __name__ == "__main__":
    rep_num = 0
    h5_path_regx = "out*/*.h5"
    track = loadTrack(rep_num, h5_path_regx)

    ener, mcc = [], []

    for i in track:
        ener.append(i[1])
        mcc.append(i[2])
        
    print pearsonr(mcc, ener)



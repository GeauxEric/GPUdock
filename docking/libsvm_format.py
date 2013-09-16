"""
#################################################################################
concat the output files and transform to the format demanded by libsvm
#################################################################################
"""
import pandas as pd
import numpy as np
import suffix

df = pd.DataFrame

def Concat(concat_lst):
    """load a list of file into df, concat them into one df"""
    all_data = []
    for file in concat_lst:
        all_data.append(pd.read_csv(file))
    data = pd.concat(all_data)

    return data

def reFormat(data):
    """reformat the data frame to the required format, return the new format"""
    data = data.sort(columns='mcc', ascending=False)
    data.index = range(len(data.index))

    true_lst = df(np.where(data['mcc'] > 0.5, -1, 1))
    true_lst.columns = ['bool']

    data = data.drop(['mcc','Unnamed: 0'], axis=1)
    old_columns = data.columns
    for num, name in enumerate(old_columns):
        tmp_df = list(data[name])
        tmp_df = df([str(num+1)+':'+("%.11f" % i) for i in tmp_df])
        tmp_df.columns=[str(num+1)]
        data = pd.concat([data, tmp_df], axis=1)
    data = data.drop(old_columns,axis=1)
    data = pd.concat([data, true_lst], axis=1)
    # rearrange the columns
    data = data.reindex(columns=pd.Index(['bool']).append(data.columns - ['bool']))

    return data

def main(rootdir, extd, output_file = 'svm_in'):
    """
    read multiple xxx-DecoyRepEner.txt files, concat them
    convert to the format that libsvm accept
    """
    files_to_read = suffix.getAbsPaths(rootdir, extd)
    data = Concat(files_to_read)
    data = reFormat(data)
    data.to_csv(output_file, sep=' ',index=False, header=False)

if __name__ == '__main__':
    import sys,os
    extd = ' '

    if len(sys.argv) == 3:
        extd = sys.argv[1]
        output_file = sys.argv[2]
        rootdir = os.getcwd()       # current dir by default
    elif len(sys.argv) == 4:
        rootdir = sys.argv[1]
        extd = sys.argv[2]
        output_file = sys.argv[3]
    else:
        print "usage:   "
        print "python libsvm_format.py <data_dir> <extd-name> <output file name>"
        sys.exit(1)

    main(rootdir, extd, output_file)

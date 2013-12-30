import sys, os
import subprocess
import glob
import pandas

df = pandas.DataFrame

# def loadDt(ifn, sep=','):
#     """load the file into DataFrame
    
#     Arguments:
#     - `ifn`: input file path
#     """
#     print "loading ...\t\t\t", ifn
#     # import ipdb; ipdb.set_trace()
#     dt = pd.read_csv(ifn)
#     # dt = pd.read_csv(ifn, header=True, sep=sep)
#     return dt

    
# def removeDup(dt):
#     """remove the duplicating rows from the table, in-place 
    
#     Arguments:
#     - `dt`: DataFrame of the data
#     """
#     print "removing duplicates ..."
#     dt = dt.drop_duplicates()
    
# def printTotal(dt, ifn):
#     """print the total energy, whose name is "total"
#     """
#     total = df(dt['total'])
#     removeDup(total)
#     ofn = ifn
#     total.to_csv('ifn')

def main():
    pattern = './output_*/a_*'
    outputs = glob.glob(pattern)
    hdf_path = outputs[0]  # in testing, only the first hdf file tested
    hdf_path = os.path.abspath(hdf_path)

    base = os.path.splitext(hdf_path)[0] 
    analysis_path = base + '.csv'

    print "pandas loading \t\t\t", analysis_path

    dt = pandas.read_csv(analysis_path)

    print "removing duplicates ..."
    dt = dt.drop_duplicates()
    
    ofn = base + '_total.csv'
    print "total enregy\t\t", ofn
    total = df(dt['total'])
    total.to_csv(ofn, index=False)


    # import ipdb; ipdb.set_trace()
    # load the plain text into pandas DataFrame
    # dt = loadDt(analysis_path)
    # dt = pd.read_csv(analysis_path)

    # # remove duplicates
    # removeDup(dt)

    # # print the total energy
    # outfn = base + '_total.csv'
    # printTotal(dt, outfn)

if __name__ == "__main__":
    main()

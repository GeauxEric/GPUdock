import sys, os
import subprocess
import glob


def readHdf(analysisProgram='./analysis', hdf_path=' ', analysis_path=' '):
    """read the Hdf log using cProgram
    
    Arguments:
    - `analysisProgram`: program to read the hdf data
    """
    print 'analysis C program\t\t', analysisProgram
    print 'hdf path\t\t\t', hdf_path
    print 'analysis output\t\t\t', analysis_path

    cmd = [analysisProgram, hdf_path]
    with open (analysis_path, 'wb') as out:
        subprocess.Popen(cmd, stdout = out)
    
    print "reading HDF file completed ..."


def main():
    # read the hdf file into plain text
    pattern = './output_*/a_*'
    outputs = glob.glob(pattern)
    hdf_path = outputs[0]  # in testing, only the first hdf file tested
    hdf_path = os.path.abspath(hdf_path)

    base = os.path.splitext(hdf_path)[0] 
    analysis_path = base + '.csv'
    
    analysisProgram = './analysis'

    readHdf(analysisProgram=analysisProgram, hdf_path=hdf_path, analysis_path=analysis_path)
    
if __name__ == "__main__":
    main()
    

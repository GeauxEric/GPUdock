"""
suffix.py:	deal with all files with the suffix .xx in a directory tree
usage:		$python surffix.py argv[1] argv[2]
extd:		the surffix: .xx
rootdir:	the directory tree
"""

import os, sys

def getPaths(rootdir, extd):
    """return a list of the absolute path of the files in rootdir end with extd"""
    fullnames = []
    for (thisDirLevel, subsHere, filesHere) in os.walk(rootdir):
        for filename in filesHere:
            if filename.endswith(extd):
                fullname = os.path.abspath(filename)
                fullnames.append(fullname)
    return fullnames

if __name__ == '__main__':
    # current dir by default
    import sys,os
    rootdir = os.getcwd() if len(sys.argv) == 2 else sys.argv[1]
    extd = ' '

    if len(sys.argv) == 2:
            extd = sys.argv[1]
    elif len(sys.argv) == 3:
            rootdir = sys.argv[1]
            extd = sys.argv[2]
    else:
            print "too few argvs..."
            sys.exit(1)		# exit FAILURE

    find = getPaths(rootdir, extd)
    print len(find),extd, 'found'

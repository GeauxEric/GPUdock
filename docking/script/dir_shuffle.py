"""
dir_shuffle.py:		iterate sub directories and return their names shuffled
usage:			$python subdir.py argv[1] 
rootdir:		the directory tree
"""

import os, sys
import random		

# current dir by default or user-specified
rootdir = os.getcwd() if len(sys.argv) == 1 else sys.argv[1]		 
# empty container
ls = []

valid = 0
for (thisDirLevel, subsHere, filesHere) in os.walk(rootdir):
	  for dirname in subsHere:
		    if os.path.isdir(dirname):		# if dirname is a dirctory
				ls.append(dirname)		# added to the list
				valid += 1

# print valid, "directories found!"

random.shuffle(ls)						# shuffle the list
for i in ls:							# print the name out
	  print i					

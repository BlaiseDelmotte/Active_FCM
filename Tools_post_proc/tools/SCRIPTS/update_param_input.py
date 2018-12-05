#!/usr/bin/python

import os,sys,time
import shutil


# Program which take new lines in input file "fparam.in"
# to update those in test_case
# Author : Blaise Delmotte
# January 2014

# Indices of lines to insert from file_orig
# SUBSTRACT -1 TO LINE NUMBERS
orig_lines_start = 184
orig_lines_end = 190

# Indices of of line to start insertion in file_dest
# SUBSTRACT -1 TO LINE NUMBERS
dest_lines_start = 184

# Defines the file where to extract data (file orig)
name_file_orig = '../../trunk/param.in'

# Defines the file where to write data (file dest)
name_file_dest = 'param.in'

# Defines the updated dest file  (file dest new)
name_file_dest_new = 'param_new.in'

# Check that this file exists
if os.path.isfile(name_file_orig) == False:
    print name_file_orig+' does not exist.'
    sys.exit(-1)
    
file_orig = open(name_file_orig,'r')

lines_orig = file_orig.readlines()

print len(lines_orig)

    
# Local variables...
dirname = "./"
folder = os.listdir(dirname)

for i in folder:
	if (os.path.isdir(i)):			
		if (i=='.svn'):
			print 'Exclude folder "',i, '"'
		else:
			print i
			print os.listdir(i)
			os.chdir(i)
			# Check that the file where to write data exists (dest file)
			if os.path.isfile(name_file_dest) == False:
				print name_file_dest,' does not exist. '
				sys.exit(-1)				
			file_dest = open(name_file_dest,'r')
			lines_dest = file_dest.readlines()
			file_dest.close()
			
			print len(lines_dest)
			
			if (len(lines_dest)<len(lines_orig)):
			
				print lines_dest
				
				k = 0
				
				for j in range(orig_lines_start,orig_lines_end):
					lines_dest.insert(dest_lines_start+k,lines_orig[j])
					k=k+1
					
				file_dest_new = open(name_file_dest_new,'w')
				file_dest_new.writelines(lines_dest)
				file_dest_new.close()
				
				os.remove(name_file_dest)
				os.rename(name_file_dest_new,name_file_dest)
			os.chdir("../")
				
			
			
			
			

#!/usr/bin/python

import os,sys,time
import shutil, glob
from subprocess import call

def message(i):
    print i
 
if __name__ == '__main__':
    print 'sys.argv: ', sys.argv
    if len(sys.argv) > 1:
        message(sys.argv[1])
    else:
        message('No argument')


To_add = int(sys.argv[1])

srcdir = './'
srcfile_bin = glob.iglob(os.path.join(srcdir, "*.bin"))



for file in srcfile_bin:
 file_split = file.split('.')
 file_split2 = file_split[1].split('t')
 time_file = file_split2[1]
 time_file_new = To_add + int(time_file)
 str_time_file_new = '0000000'+ str(time_file_new)
 length_num = len(str_time_file_new)
 number_new = str_time_file_new[length_num-8:length_num]
 #print time_file
 file_split3 = file_split2[0].split('/')
 new_name = file_split3[1] + 't' + number_new + '.bin'
 #print new_name
 os.rename(file,new_name)
 #raw_input()

 

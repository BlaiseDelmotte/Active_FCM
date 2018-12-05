#!/usr/bin/python

import os,sys,time
import shutil, glob
from subprocess import call

srcdir = './'
srcfile_ascii = glob.iglob(os.path.join(srcdir, "Frame.*.png"))


for file in srcfile_ascii:
 print file
 file_split = file.split('.')
 file_split[2] = '000000'+ file_split[2]
 length_num = len(file_split[2])
 number_new = file_split[2][length_num-6:length_num]
 new_filename = './Frame_' + number_new + '.png'
 print new_filename
 os.rename(file,new_filename)
 
 

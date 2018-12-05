#!/usr/bin/python

import os,sys,time
import shutil, glob
from subprocess import call
#~ import subprocess as sub

def message(i):
    print i


if __name__ == '__main__':
    print 'sys.argv: ', sys.argv
    if len(sys.argv) > 1:
        message(sys.argv)
    else:
        message('No argument: voronoi for all files')


# Program which compute Voronoi Tesselation
# for given ascii files 
# Author : Blaise Delmotte
# April 2014



###### REQUIRES REFLEXION : choose only files between save_start
###### and save_end
if len(sys.argv) > 1:
	# Defines the starting save
	save_start = sys.argv[1]
	# Defines the ending save
	save_end = sys.argv[2]

	message(save_start)
	message(save_end)


# Defines the file where to extract data (file orig)
list_file = []
srcdir = './'
if len(sys.argv) == 1:
	srcfile_ascii = glob.iglob(os.path.join(srcdir, "FCM_PART_POS_VORO_t*.dat"))
	for file in srcfile_ascii:
		list_file.append(file)

else:
	srcfile_ascii = glob.iglob(os.path.join(srcdir, "FCM_PART_POS_VORO_t*.dat"))
	for file in srcfile_ascii:
		file_split = file.split('0',1)
		file_split = file_split[1].split('.',1)
		file_split = file_split[0]
		
		if int(file_split)>=int(save_start) and int(file_split)<=int(save_end):
			list_file.append(file)


# Read fcm_run.info to get domain size
name_file_info = 'fcm_run.info'
file_info = open(name_file_info,'r')
lines_info = file_info.readlines()

domain_size = float(lines_info[0])
rad_part = float(lines_info[14])
max_size = lines_info[16]
max_size = max_size.split('0 ',1)
min_size = float(max_size[1])
max_size = float(max_size[0])


size_adim = domain_size/rad_part*max_size/min_size

# Calls Voro++ with arguments:
#    - "-g" : gnuplot files to see tesselation,
#    - "-c" "%i %x %y %z %s %v" : custom format saving with ID / POS / #FACES /VOLUME
#    - "-p" : periodic boundary conditions
for file in list_file:
	cmd = 'voro++ -g -c "%i %x %y %z %s %v" -p 0 ' + str(size_adim) + ' 0 ' + str(size_adim) + ' 0 ' + str(size_adim) + ' ' + file
	print cmd
	os.system(cmd) 
	file_dat_vol = file + '.vol'
	file_voro = file_dat_vol.replace('.dat.vol','.voro')
	os.rename(file_dat_vol,file_voro)
	file_dat_gnu = file + '.gnu'
	file_gnu = file_dat_gnu.replace('.dat.gnu','.gnu')
	os.rename(file_dat_gnu,file_gnu)
	

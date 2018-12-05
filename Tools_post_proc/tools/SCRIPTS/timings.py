import os,sys,time
import tarfile
import shutil, glob
import subprocess

srcdir = './'
Rfolders = glob.iglob(os.path.join(srcdir, "R*"))
name_file_dest = "ns3d_20cores.log"
num_line = 615

k=0
timing = []
for folders in Rfolders:
	if os.path.isfile(folders) == False:
		k = k+1
		os.chdir(folders)
		file_dest = open(name_file_dest,'r')
		lines_dest = file_dest.readlines()
		line_final_time = lines_dest[num_line-1]
		line_final_time = line_final_time.split(':')
		line_final_time = line_final_time[1].split('s')
		timing.append(line_final_time[0])
		os.chdir('../')

f= open('List_timings.txt','w')
for i in range(0,len(timing)-1):
	f.write(timing[i]+"\n")
f.close()

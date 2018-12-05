#!/usr/bin/python

import os,sys,time
import tarfile
import shutil, glob
import subprocess

srcdir = './'
cmd_connec = "./treat_part"
os.system(cmd_connec)
files_to_treat = glob.iglob(os.path.join(srcdir, "CONNECTIVITY_t*"))

for file in files_to_treat :
	file_split = file.split('t')
	file_split = file_split[1].split('.')
	time = file_split[0]
	output_file = "2CLIQUES_t"+time+".dat"
	cmd_cliques = "python PercAlgo.py " + file + " 2 " + output_file 
	print cmd_cliques
	os.system(cmd_cliques)
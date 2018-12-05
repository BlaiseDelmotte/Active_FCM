#!/usr/bin/python

import os,sys,time
import tarfile
import shutil, glob
from subprocess import call

srcdir = './'
srcfile_tar = glob.iglob(os.path.join(srcdir, "results.*.tar"))

k=0
for file in srcfile_tar:
 k = k+1
 os.makedirs("R"+str(k))
 os.makedirs("R"+str(k)+"/Data_files/")
 os.makedirs("R"+str(k)+"/Ascii_files/")
 file_split = file.split('.tar')
 file_split = file_split[0].split('./')
 filename = file_split[1]
 tar = tarfile.open(file)
 tar.extractall()
 tar.close()
 folder_to_move = "./workgpfs/rech/uxf/ruxf002/"+filename
 folder_dest = "R"+str(k)+"/Data_files/"
 files_to_move = glob.iglob(os.path.join(folder_to_move, "*.*"))
 print folder_to_move
 print folder_dest
 for file_to_move in files_to_move:
  if (file_to_move=="."):
   print 'Exclude folder "',file_to_move, '"'
  else:
   shutil.move(file_to_move,folder_dest)
   
 shutil.rmtree("./workgpfs")
 
 shutil.copy("./mpi2ascii",folder_dest)
 shutil.copy("./convert_ascii_choice.input",folder_dest)
 
 os.chdir(folder_dest)
 os.system("./mpi2ascii")
 folder_ascii_dest = ("../Ascii_files")
 ascii_files = glob.iglob(os.path.join("./", "*.dat"))
 
 for ascii_to_move in ascii_files:
  shutil.move(ascii_to_move,folder_ascii_dest)
  
 in_files = glob.iglob(os.path.join("./", "*.in*"))
 for in_to_move in in_files:
  shutil.copy(in_to_move,folder_ascii_dest)
  
 os.chdir("../../")
 

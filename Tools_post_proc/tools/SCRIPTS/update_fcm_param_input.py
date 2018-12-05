#!/usr/bin/python

import os,sys,time
import shutil


# Program which take new lines in input file "fcm_param.in"
# to update those in test_case
# Author : Blaise Delmotte
# June 2013

# Indices of lines to insert from file_orig
orig_lines_start = 67
orig_lines_end = 68

# Indices of of line to start insertion in file_dest
dest_lines_start = 67

# Defines the file where to extract data (file orig)
name_file_orig = '../../trunk/fcm_param.in'

# Defines the file where to write data (file dest)
name_file_dest = 'fcm_param.in'

# Defines the updated dest file  (file dest new)
name_file_dest_new = 'fcm_param_new.in'

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
			
			#~ if (len(lines_dest)<len(lines_orig)):
			
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
				
			
			
			
			
			
			
			
			
		
		
		
		
	

#~ 
#~ # Check that the file where to write data exists (dest file)
#~ if os.path.isfile('input.in') == False:
    #~ print 'input.in does not exist. exiting.'
    #~ sys.exit(-1)
#~ # Read the source file
#~ file = open(cfile,'r')
#~ lines = file.readlines()
#~ file.close()
#~ 
#~ # Number of params
#~ number_params=len(lines[0].split(','))
#~ 
#~ # Number of n-uplet for params
#~ number_nuplet=len(lines)
#~ 
#~ # Lines to modify in the dest file
#~ name_params=[]
#~ name_params=lines[0].split(',')
#~ 
#~ for k in range(number_params):
	#~ name_params[k]=str(name_params[k]).replace('\n','')
#~ 
#~ 
#~ print name_params
#~ # Once first line treated, delete it
#~ lines.pop(0)
#~ 
#~ # Lines to modify in the dest file
#~ lines_input=[]
#~ lines_input=lines[0].split(',')
#~ 
#~ for k in range(number_params):
	#~ lines_input[k]=int(str(lines_input[k]).replace('\n',''))-1
	#~ 
#~ # Once first line treated, delete it
#~ lines.pop(0)
#~ 
#~ print lines
#~ 
# Number of lines in histo_nit.txt
#~ #num_lines_histo=25
#~ 
#~ 
#~ 
#~ for i in range(number_nuplet):
	#~ 
	#~ print i
#~ 
    #~ # Save the previous version of input.dat
	#~ os.system('cp input.in input.in.prev')
#~ 
    #~ # Read then modify the lines of input.in
	#~ file = open('input.in.ref','r')
	#~ lines_dest = file.readlines()
	#~ file.close()
    #~ 
    #~ # Extract the new parameter values
	#~ values = lines[i].split(',')
    #~ 
    #~ # Copy the new parameter values
	#~ value_str = []
    #~ 
	#~ for k in range(number_params):
		#~ value_str.append( str(values[k]).replace('\n',''))
		#~ 
		#~ # before EOL, need to add "\n"
		#~ lines_dest[lines_input[k]] = value_str[k]+"\n"
	#~ 
#~ 
	#~ # Write the new version of input.dat
	#~ file = open('input.in','w')
	#~ for j in range(len(lines_dest)):
		#~ file.write(lines_dest[j])
	#~ 
	#~ file.close()
	#~ 
	#~ # Compilation + execution
	#~ print '\n Executing input.in with params = '
	#~ print value_str
	#~ 
	#~ 
	 #~ 
	#~ 
	#~ # Create the associate folder to the n-uplet's value
	#~ folder_name = ""
    #~ 
	#~ for k in range(number_params):
		#~ if k==0:
			#~ folder_name = folder_name + name_params[k] + '_' + value_str[k].replace('.','_')
		#~ else:
			#~ folder_name = folder_name + '_' + name_params[k] + '_' + value_str[k].replace('.','_')
	#~ 
	#~ os.makedirs(folder_name)
	#~ 
	#~ # Copy the new 'input.in' in this folder
	#~ src = "input.in"
	#~ dst = folder_name+"/"+src
	#~ shutil.copyfile(src, dst)
	#~ 
	#~ os.chdir(folder_name)
	#~ os.system("pwd")
	#~ os.system("/bin/csh -i -c 'run4_SpSt_WorkDev'")
	#~ 
	#~ os.chdir("../")
	
	#~ os.system('source .bashrc; shopt -s expand_aliases; run4_SpSt_WorkDev')
	
	#~ 
	#~ run4.communicate()


    
    
    
    
		
    

    
    
    
    
#~ 
    #~ # Save the last version of input.dat
    #~ os.system('cp input.dat input.dat.last')
    #~ 
    #~ # Write Histo_nit.txt files for each external loop
    #~ for k in range(n_ext):
	#~ src = "Histo_"+str(k+1).zfill(2)+".txt"
	#~ file = open(src,'r')
   	#~ lines_histo = file.readlines()
   	#~ file.close()
    	#~ dst = folderName+"/"+src
	#~ shutil.copyfile(src, dst)
	#~ if len(lines_histo)<num_lines_histo:
		#~ break
	#~ elif float(lines_histo[line_dist])<eps:
		#~ break
#~ 
      #~ 
      #~ 
#~ 
#~ 

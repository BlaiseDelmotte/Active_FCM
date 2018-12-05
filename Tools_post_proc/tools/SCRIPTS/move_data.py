import os,sys,time
import shutil, glob
 
def message(i):
    print i
 
if __name__ == '__main__':
    print 'sys.argv: ', sys.argv
    if len(sys.argv) > 1:
        message(sys.argv[1])
    else:
        message('No argument')
      

# Defines the folder where to copy data (folder dest)
folder_dest = sys.argv[1]
if os.path.isdir(folder_dest) == False:
 print folder_dest+' does not exist.'
 print 'I create it.'
 os.makedirs(folder_dest)

# Defines the file where to extract data (file orig)
srcdir = './'
srcfile_in = glob.iglob(os.path.join(srcdir, "*.in"))
srcfile_info = glob.iglob(os.path.join(srcdir, "*.info"))
srcfile_bin = glob.iglob(os.path.join(srcdir, "*.bin"))
srcfile_end = glob.iglob(os.path.join(srcdir, "*.end"))


# Check that these file exist and copy them
for file in srcfile_in:
    if os.path.isfile(file):
        shutil.copy2(file, folder_dest)

for file in srcfile_info:
    if os.path.isfile(file):
        shutil.copy2(file, folder_dest)
        
for file in srcfile_bin:
    if os.path.isfile(file):
        shutil.copy2(file, folder_dest)
        
for file in srcfile_end:
    if os.path.isfile(file):
        shutil.copy2(file, folder_dest)



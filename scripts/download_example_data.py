# This example script downloads the ALFALFA grid 1044+13 from Zenodo.
# These data are packaged exactly as if downloaded from NRAO.
# The script will create a "data" directory in the main directory,
# then download the tar file, and finally untar it.

import os, subprocess, urllib

#Check which is the current directory
cwd = os.getcwd()+'/'
if ('scripts' in cwd[-8:-1]) or ('notebooks' in cwd[-10:-1]) :
    directory_path = cwd + '../data'
elif 'AALegacy' in cwd[-9:-1]:
    directory_path = cwd + '/data'
else:
    print("This script does not appear to be running in the scripts directory.")
    print("It will now exit without completing.")
    print("Please re-run in the scripts directory.")
    exit(-1)

#Make new directory
if not os.path.exists(directory_path):
    os.makedirs(directory_path)
    print(f"Directory '{directory_path}' created successfully!")
else:
    print(f"Directory '{directory_path}' already exists.")

#Downlaod file
#######

#Untar the file
run_result = subprocess.run(['tar','-xvf',
                             directory_path+'/NRAO_archive_A2010_20250627-074204.tar',
                             '-C',directory_path+'/'],capture_output=True, text=True)
print(run_result.stdout)
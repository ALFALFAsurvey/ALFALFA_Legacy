# This example script downloads the ALFALFA grid 1044+13 from Zenodo.
# These data are packaged exactly as if downloaded from NRAO.
# The script will create a "data" directory in the main directory,
# then download the tar file, and finally untar it.

import os, subprocess, urllib.request

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
    exit()

#Make new directory
if not os.path.exists(directory_path):
    os.makedirs(directory_path)
    print(f"Directory '{directory_path}' created successfully!")
else:
    print(f"Directory '{directory_path}' already exists.")

#Downlaod file
def download_progress_hook(count, block_size, total_size):
    if total_size > 0:
        percentage = min(100, int(100 * count * block_size / total_size))
        print(f"Downloading: {percentage}%", end='\r')
    else:
        print(f"Downloading: {count * block_size / (1024*1024):.2f} MB downloaded", end='\r')
        
def download_file(url, filename):
    try:
        print(f"Starting download of {url} to {filename}...")
        urllib.request.urlretrieve(url, filename, reporthook=download_progress_hook)
        print("\nDownload complete!")
    except Exception as e:
        print(f"\nError during download: {e}")
        if os.path.exists(filename):
            print("File already exists.")
            print("Please (re)move it and re-run the script.")
            exit()
            
download_file("https://zenodo.org/records/16592112/files/NRAO_archive_A2010_20250730-112528.tar",
              directory_path+'/NRAO_archive_A2010_20250730-112528.tar')

#Untar the file
run_result = subprocess.run(['tar','-xvf',
                             directory_path+'/NRAO_archive_A2010_20250730-112528.tar',
                             '-C',directory_path+'/'],capture_output=True, text=True)
print(run_result.stdout)

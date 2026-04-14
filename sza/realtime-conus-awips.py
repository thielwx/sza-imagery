#============================================
# Code that checks for CMI data, and if avaialble it produces SZA imagery as netCDF files
# Currently this is designed to run for only CONUS data (East/West), but it can be adapted in the future
# to run for MESO data too.
#
# Created: April 2026
# Author: Kevin Thiel (kevin.thiel@ou.edu)
#===========================================

import sza_calc as sza #'Library' of functions I pull from
from datetime import datetime, timedelta, UTC
import os
import sys
import yaml
import multiprocess as mp
from glob import glob
import itertools


# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Change the current working directory to the script's directory
os.chdir(script_dir)

#Getting the current time
start_datetime = datetime.now(UTC)
#start_datetime = datetime(2025,5,29,0,37,) #DEVMODE
cur_file_time = sza.file_time_converter(start_datetime-timedelta(minutes=5))
y, m, d, doy, hr, mi = sza.datetime_converter(cur_file_time)
dload_time = y+'-'+m+'-'+d+'T'+hr+':'+mi+':00Z'
file_time = 's'+y+doy+hr+mi

#start_datetime = datetime.now(UTC) #DEVMODE

#Importing the appropriate yaml file
with open('realtime-conus-awips.yaml', 'r') as f:
    sfile = yaml.safe_load(f)

#Loading the universal variables
channels = sfile['channels']
slots = sfile['slots']
scenes = sfile['scenes']

input_file_loc = sfile['input-file-loc']
output_loc = sfile['output-loc']

timeout_limit = sfile['timeout-limit']
cpu_max = sfile['cpu-max']



#===========================
# Phase 1: Processing the data for AWIPS with SZA channels
#===========================

#Looking for the input files that we'll process
input_file_list = []
for combo in itertools.product(scenes, channels, slots):
    scene = combo[0]
    channel = combo[1]
    slot = combo[2]
    
    search_str = input_file_loc + 'OR_ABI-L2-CMIP' + scene + '-M?' + channel + '_' + slot + '_' + file_time + '*.nc'

    file_search = glob(search_str)

    #First check that there is data
    if len(file_search)==1:
        
        check_file = file_search[0]
        
        #Starting second check that we're not constantly overwriting data
        idx = check_file.find('CMIP') #Finding the CMIP tag
        sza_file_search = check_file[:idx]+'SZA'+check_file[idx:]
        out_file_search = glob(sza_file_search) #inserting SZA name to hunt for matching file
        
        #If there's an input file, but not already an output file, put the file in the list to be processed
        if len(out_file_search)==0:
            input_file_list.append(file_search[0])

        
#If there's no files to process, we're done!
if len(input_file_list)==0:
    sys.exit(0)

#Creating lists to put into the multiprocessing
output_list = [output_loc for i in range(len(input_file_list))]

if len(input_file_list)<cpu_max:
    n_cpu = len(input_file_list)
else:
    n_cpu = cpu_max


#Driver funciton that processes the data in parallel
def realtime_driver(temp_file, output_loc):
    sza.sza_io(temp_file, output_loc)


#sending the file names to be processed
if __name__ == "__main__":
    with mp.Pool(processes=n_cpu) as p:
        async_result = p.starmap_async(realtime_driver,zip(input_file_list, output_list))

        #Ensuring we don't get a stuck process
        try:
            results = async_result.get(timeout=timeout_limit)
        except mp.TimeoutError:
            print ('Timeout error: '+dload_time)

        p.close()
        p.join()


print('**************'+str(start_datetime)+'**************')
print ('Finshed, runtime: '+str(datetime.now(UTC)-start_datetime))

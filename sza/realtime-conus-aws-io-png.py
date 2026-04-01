#============================================
# Code that downloads CMI data from AWS, then produces sza output files
#
# Created: March 2026
# Author: Kevin Thiel (kevin.thiel@ou.edu)
#===========================================

import sza_calc as sza
from datetime import datetime, timedelta, UTC
import os
import yaml
import multiprocess as mp
import fnmatch
import netCDF4 as nc

# Import the locator and datasource according to your desired product
from goesdl.goesr import GOESProductLocatorABIPP
from goesdl.datasource import DatasourceAWS
from goesdl.downloader import Downloader
  

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Change the current working directory to the script's directory
os.chdir(script_dir)

#Getting the current time
start_datetime = datetime.now(UTC)
cur_file_time = sza.file_time_converter(start_datetime-timedelta(minutes=5))
y, m, d, doy, hr, mi = sza.datetime_converter(cur_file_time)
dload_time = y+'-'+m+'-'+d+'T'+hr+':'+mi+':00Z'
file_time = 's'+y+doy+hr+mi

print('**************'+str(start_datetime)+'**************')

#Importing the appropriate yaml file
with open('realtime-conus-aws-io.yaml', 'r') as f:
    sfile = yaml.safe_load(f)

#Loading the universal variables
channels = sfile['channels']
slots = sfile['slots']
temp_file_loc = sfile['temp-file-loc']
output_loc = sfile['output-loc']
sza_threshold = sfile['sza-threshold']
timeout_limit = sfile['timeout-limit']

rgb_extra_channels = sfile['rgb-channels']
fig_output_loc = sfile['fig-output-loc']


#===========================
# Phase 1: Processing the data for AWIPS with SZA channels
#===========================

#Downloading the data that needs the solar zenith adjustment
temp_file_list  = sza.aws_downloader(slots, channels, temp_file_loc, dload_time, dload_time)

#Creating lists to put into the multiprocessing
output_list = [output_loc for i in range(len(temp_file_list))]
szath_list = [sza_threshold for i in range(len(temp_file_list))]


#Driver funciton that processes the data in parallel
def realtime_driver(temp_file, output_loc, sza_threshold):
    sza.sza_io(temp_file, output_loc, sza_threshold)


#sending the file names to be processed
if __name__ == "__main__":
    with mp.Pool(processes=12) as p:
        async_result = p.starmap_async(realtime_driver,zip(temp_file_list, output_list, szath_list))

        #Ensuring we don't get a stuck process
        try:
            results = async_result.get(timeout=timeout_limit)
        except mp.TimeoutError:
            print ('Timeout error: '+dload_time)

        p.close()
        p.join()


#================================
# Phase 2: Making the images for CIRA slider
#================================

#Downlaoding the extra RGB channel(s) from AWS
#temp_file_rgb_list  = sza.aws_downloader(slots, rgb_extra_channels, temp_file_loc, dload_time, dload_time)


#r_files = []
g_files = []
#b_files = []

realtime_files = os.listdir(output_loc)

#Loops for collecting the rgb files
for slot in slots:
    #Red 
    #pattern = '*C13_'+slot+'_'+file_time+'*.nc'
    #r_files.extend([item for item in temp_file_rgb_list if fnmatch.fnmatch(item, pattern)])
    #Green
    pattern = '*C02_'+slot+'_'+file_time+'*.nc'
    g_files.extend([output_loc+item for item in realtime_files if fnmatch.fnmatch(item, pattern)])
    #Blue
    #pattern = '*C05_'+slot+'_'+file_time+'*.nc'
    #b_files.extend([output_loc+item for item in realtime_files if fnmatch.fnmatch(item, pattern)])


#Creating lists to put into the multiprocessing
fig_output_loc_list = [fig_output_loc for i in range(len(slots))]
file_time_list = [file_time for i in range(len(slots))]

#def realtime_image_driver(red, green, blue, output_loc, slot, time_string):
def realtime_image_driver(green, output_loc, slot, time_string):
    #Plotting the ch2 data first
    sza.realtime_slider_ch02(green, slot, output_loc, time_string)

    #NO SZA ADJUSTMENT SINCE WERE PULLING FROM THE SZA FILES THEMSELVES
    #rgb = sza.rgb_creator_dcpd(nc.Dataset(red), nc.Dataset(green), nc.Dataset(blue), scene='CONUS', sza_adjustment=False) 

    #Plotting the rgbs
    #sza.rgb_plotter_v2(rgb, slot, output_loc, time_string)


#sending the file names to be processed
if __name__ == "__main__":
    with mp.Pool(processes=12) as p:
        #async_result = p.starmap_async(realtime_image_driver,zip(r_files, g_files, b_files, fig_output_loc_list, slots, file_time_list))
        async_result = p.starmap_async(realtime_image_driver,zip(g_files, fig_output_loc_list, slots, file_time_list))

        #Ensuring we don't get a stuck process
        try:
            results = async_result.get(timeout=timeout_limit)
        except mp.TimeoutError:
            print ('Timeout error: '+dload_time)

        p.close()
        p.join()


#================================
# Wrapping up the script
#================================

#The final step: Removing the temporary files
for temp_file in temp_file_list:
    os.remove(temp_file)
#for temp_file in temp_file_rgb_list:
#    os.remove(temp_file)

print ('Finshed, runtime: '+str(datetime.now(UTC)-start_datetime))

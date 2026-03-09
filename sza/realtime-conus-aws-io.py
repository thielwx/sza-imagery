#============================================
# Code that downloads CMI data from AWS, then produces sza output files
#
# Created: March 2026
# Author: Kevin Thiel (kevin.thiel@ou.edu)
#===========================================

import sza_calc as sza
from datetime import datetime, timedelta
import os
import yaml
import multiprocess as mp

# Import the locator and datasource according to your desired product
from goesdl.goesr import GOESProductLocatorABIPP
from goesdl.datasource import DatasourceAWS
from goesdl.downloader import Downloader


# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Change the current working directory to the script's directory
os.chdir(script_dir)

#Getting the current time
start_datetime = datetime.now()
cur_file_time = file_time_converter()
y, m, d, doy, hr, mi = sza.datetime_converter(cur_file_time)
dload_time = y+'-'+m+'-'+d+'T'+hr+':'+mi+':00Z'


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

#Empty file list we'll use to know where our temporary files are
temp_file_list = []

#Looping in case we're pulling from multiple satellites (GOES-East/-West)
for slot in slots:
    # Initialize the product locator for GOES-R Series (set your desired product)
    locator = GOESProductLocatorABIPP("CMIP", "C", channels, slot)
    datasource = DatasourceAWS(locator, cache=600)

    # Initialize the downloader with the locator and datasource
    downloader = Downloader(
        datasource=datasource,
        locator=locator,
        repository=temp_file_loc,
    )

    #Downloading the files from AWS at the time specified
    t_files = downloader.download_files(start=dload_time, end=dload_time)
    
    #Adding the files to a list that we'll pull from
    t_files_new = [temp_file_loc+f for f in t_files]
    temp_file_list = temp_file_list.extend(t_files_new)


#sending the file names to be processed
if __name__ == "__main__":
    with mp.Pool(processes=12) as p:
        p.starmap_async(realtime_driver,temp_file_list)

        #Ensuring we don't get a stuck process
        try:
            results = async_result.get(timeout=timeout_limit)
        except mp.TimeoutError:
            print ('Timeout error: '+dload_time)

        p.close()
        p.join()



#Driver funciton that processes the data in parallel
def realtime_driver(temp_file):
    global output_loc
    global sza_threshold
    sza.sza_io(temp_file, output_loc, sza_threshold)

    os.remove(temp_file)


    












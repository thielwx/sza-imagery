#============================================
# Code that removes old data (2+ hrs) when running realtime-conus-aws-io.py
# This script should be run hourly
#
# Created: March 2026
# Author: Kevin Thiel (kevin.thiel@ou.edu)
#===========================================

import sza_calc as sza
from datetime import datetime, timedelta
import os
import shutil
from glob import glob
import yaml


# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Change the current working directory to the script's directory
os.chdir(script_dir)

#Importing the appropriate yaml file and needed constants
with open('realtime-conus-aws-io.yaml', 'r') as f:
    sfile = yaml.safe_load(f)

temp_file_loc = sfile['temp-file-loc']
output_loc = sfile['output-loc']

#Getting the current time
start_datetime = datetime.now()

tdelta_hr = timedelta(hours=3)
tdelta_1day = timedelta(hours=24)

#Only run the daily wiper when more than three hours into a day
if int(datetime.strftime(start_datetime,'%H'))>=3:
    #Getting the previous day string needed to check    
    t = start_datetime - tdelta_1day
    y, m, d, doy, hr, mi = sza.datetime_converter(t)
    day_str_temp = temp_file_loc + y + '/' + doy + '/'
    
    #If the path still exists, remove yesterday's directory
    if os.path.exists(day_str_temp):
        try:
            shutil.rmtree(day_str_temp)        
        except OSError as e:
            print(f"Error removing directory: {e}")


#Getting the previous hours string that we can pull data from
t = start_datetime - tdelta_hr
y, m, d, doy, hr, mi = sza.datetime_converter(t)
hr_str_temp = output_file_loc + '*s' + y + doy + hr

file_list = glob(hr_str_temp)

#If there's files, delete them
if len(file_list)>0:
    for f in file_list:
        try:
            os.remove(f)
        except FileNotFoundError:
            print(f"File '{file_path}' not found.")
        except PermissionError:
            print(f"Permission denied to delete the file '{file_path}'.")
        except Exception as e:
            print (f"Error occurred while deleting the file '{f}': {e}")


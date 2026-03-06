#==============================================
# A script that's designed to run the SZA imagery in realtime for CONUS data
# The current input source is the NASA SPoRT web server
#
# Author: Kevin Thiel (kevin.thiel@ou.edu)
# Created: March 2026
#==============================================

from datetime import datetime, timedelta
import sys
import os
import sza_calc as sza
import yaml
import subprocess as sp
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import multiprocess as mp



# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Change the current working directory to the script's directory
os.chdir(script_dir)

#Getting the current time
start_datetime = datetime.now()
cur_file_time = file_time_converter()
y, m, d, doy, hr, mi = sza.datetime_converter(cur_file_time)
time_str = 's'+y+doy+hr+mi

#Importing the appropriate yaml file
with open('sza-realtime-conus.yaml', 'r') as f:
    sfile = yaml.safe_load(f)

#Loading the universal variables
channels = sfile['channels']
slots = sfile['slots']
timeout_limit = sfile['timeout-limit']

#Getting the combinations of slots and channels to pull from
combos = [(c,s) for c in channels for s in slots]



if __name__ == "__main__":
    with mp.Pool(processes=12) as p:
        p.starmap(realtime_driver,combos)
        p.close()
        p.join()




def file_time_converter(time):
    #Getting the time difference from the 0 and 5 ones place
    dt_int = int(time.strftime('%M'))%5

    if dt_int==0:
        dt_int=5
    
    #Getting the target file times from the most recent file 
    new_time = time - timedelta(minutes=int(dt_int)-1)

    return new_time

def check_file_exists_urllib(url):
    """
    Checks if a file exists on a web server using urllib.request.
    """
    # Create a request object and set the method to HEAD
    req = Request(url, method='HEAD')
    try:
        # Try to open the URL
        with urlopen(req) as response:
            #print(f"File exists. Status code: {response.getcode()}")
            return True
    except HTTPError as e:
        # Handle specific HTTP errors (e.g., 404 Not Found, 403 Forbidden)
        #print(f"File not found or access denied. Status code: {e.code}")
        return False
    except URLError as e:
        # Handle general URL errors (e.g., connection issues)
        #print(f"Connection error: {e.reason}")
        return False


def realtime_driver(channel, slot):
    global sfile
    sza_threshold = sfile['sza-threshold']
    temp_file_loc = sfile['temp-file-loc']
    output_file_loc = sfile['output-file-loc']
    orbital_slot = sfile[slot]['orbital-slot']
    sport_url = sfile[slot]['sport-url']

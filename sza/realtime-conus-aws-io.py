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


import logging
import sys


def go():
    

    return 0


# Configure logging to write to a specific file with a timestamp and level
logging.basicConfig(
    filename='/localdata/sza-realtime/realtime_io.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def main():
    try:
        now = datetime.now()
        logging.info(f"Script started at: {now}")
        # Your main script logic here

        # Get the directory where the script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # Change the current working directory to the script's directory
        os.chdir(script_dir)

        #Getting the current time
        start_datetime = datetime.now()
        cur_file_time = sza.file_time_converter(start_datetime)
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
            temp_file_list.extend(t_files_new)



        #Driver funciton that processes the data in parallel
        def realtime_driver(temp_file):
            global output_loc
            global sza_threshold
            sza.sza_io(temp_file, output_loc, sza_threshold)

            os.remove(temp_file)


        #sending the file names to be processed
        if __name__ == "__main__":
            with mp.Pool(processes=12) as p:
                async_result = p.starmap_async(realtime_driver,zip(temp_file_list))

                #Ensuring we don't get a stuck process
                try:
                    results = async_result.get(timeout=timeout_limit)
                except mp.TimeoutError:
                    print ('Timeout error: '+dload_time)

                p.close()
                p.join()

        ex_time = datetime.now()-now
        logging.info(f"Task completed successfully. {ex_time}")
    except Exception as e:
        logging.error(f"An error occurred: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()



    












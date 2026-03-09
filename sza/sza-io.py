#==============================================
# A script that takes in a netCDF file and outputs its SZA value in a new file
#
# Author: Kevin Thiel (kevin.thiel@ou.edu)
# Created: March 2026
#==============================================

import sza_calc as sza
import os
import sys
import netCDF4 as nc
import shutil
import numpy as np

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Change the current working directory to the script's directory
os.chdir(script_dir)

#Reading in input variables from user
args = sys.argv
#args = ['','../test/20260224-1411-combo/OR_ABI-L2-CMIPC-M6C02_G19_s20260551411177_e20260551413550_c20260551414061.nc', '../test/20260224-1411-combo/', '88.85'] #Devmode
input_file_str = args[1]
output_file_loc = args[2]
sza_threshold = float(args[3])

#Running in sza_calc now as a function to streamline realtime processing
sza.sza_io(input_file_str, output_file_loc, sza_threshold)

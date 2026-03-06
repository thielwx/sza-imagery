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
sza_threshold = args[3]

if not os.path.exists(output_file_loc):
    os.makedirs(output_file_loc)

#Creating the new file str to (almost) match the old one
target_str = '_ABI-L2-'
target_str_loc = input_file_str.find(target_str)
output_file_str = input_file_str[target_str_loc-2:]

#Modifying the new output file string to include SZA
target_str = 'CMIP'
target_str_loc = output_file_str.find(target_str)
output_file_str = output_file_str[:target_str_loc]+'SZA'+output_file_str[target_str_loc:]

if os.path.exists(output_file_loc+output_file_str):
    print ('WARNING: FILE ALREADY EXISTS')
    print (output_file_loc+output_file_str)
    #sys.exit()

#Creating a copy of the image file, then we'll overwrite the copy
shutil.copy(input_file_str, output_file_loc+output_file_str)

#Opening the file in write mode
dset = nc.Dataset(output_file_loc+output_file_str, 'r+')

#Getting the sza cmi values (currently using default sza_angle_threshold)
cmi_sza = sza.sza_calculator_v2_exact(dset, sza_threshold=sza_threshold)

#Overwriting the cmi data with the sza cmi values
dset.variables['CMI'][:] = cmi_sza

dset.dataset_name = output_file_str
dset.title =  dset.title+' (Solar Zenith Angle-Adjusted)'
dset.summary = dset.summary+' Original imagery has been adjusted by its solar zenith angle from https://github.com/thielwx/sza-imagery/.'

#Closing the dataset
dset.close()

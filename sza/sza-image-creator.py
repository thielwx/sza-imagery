#==============================================
# Code for displaying SZA imagery (single channel)
#
# Author: Kevin Thiel (kevin.thiel@ou.edu)
# Created: February 2026
#==============================================

import yaml
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys
import numpy as np
from datetime import datetime,timedelta
from glob import glob
import pandas as pd
import netCDF4 as nc
import sza_calc as sza
import os

#Importing the appropriate yaml file
with open('sza-image-creator.yaml', 'r') as f:
    sfile = yaml.safe_load(f)

args = sys.argv
#args = ['','20250528-tx-supercell', '1'] #Devmode
case = args[1]

#Reading in data from yaml file
input_loc = sfile[case]['input-loc']
output_loc = sfile[case]['output-loc']

stime_str = sfile[case]['start-time']
etime_str = sfile[case]['end-time']

pltname = sfile[case]['plot-name']
pbounds = sfile[case]['plot-bounds']

dt_str = sfile[case]['delta-time']

sza_switch = sfile[case]['sza-calc']
ch_num = sfile[case]['ch-num']

#Getting the time variables established
start_time = datetime.strptime(stime_str, '%Y%m%d%H%M')
end_time = datetime.strptime(etime_str, '%Y%m%d%H%M')
time_list = pd.date_range(start=start_time, end=end_time, freq=dt_str+'min')



#===============================
# FUNCTION LAND
#===============================

def datetime_converter(time):
    '''
    This function takes in a datetime object and returns strings of time features
    PARAMS:
        time: input time (datetime object)
    Returns
        y: year (4 digit number as str)
        doy: day of year (3 digit number as str)
        hr: hour of day (2 digit number as str)
        mi: minute of hour (2 digit number as str)
        d: day of month (2 digit number as str)
        m: month of year (2 digit number as str)
    '''
    y = datetime.strftime(time,'%Y') #Year
    doy = datetime.strftime(time,'%j') #Day of year
    hr = datetime.strftime(time,'%H') #Hour
    mi = datetime.strftime(time,'%M') #Minute
    d = datetime.strftime(time,'%d') #Day of month
    m = datetime.strftime(time,'%m') #Month
    
    return y, m, d, doy, hr, mi


def cmi_file_finder(input_loc,input_file_time_str):
    file_str = input_loc+'$'+input_file_time_str+'$.nc'
    file_search = glob(file_str)

    if len(file_search)==1:
        dset = nc.Dataset(file_search[0])
        data_check = True
    else:
        print ('File Seach Error')
        print (file_str)
        dset = ['NULL']
        data_check = False

    return dset, data_check


def plot_data(dset):
    geo_crs = ccrs.Geostationary(central_longitude=dset.variables['goes_imager_projection'].longitude_of_projection_origin, satellite_height=dset.variables['goes_imager_projection'].perspective_point_height, sweep_axis='x')
    extent = np.concat((dset.variables['x_image_bounds'][:],dset.variables['y_image_bounds'][::-1]), axis=0) * dset.variables['goes_imager_projection'].perspective_point_height
    return geo_crs, extent


#=====================================
# PLOT LAND
#=====================================

#Looping through each datetime
for t in time_list:

    #Arranging the time input
    y, m, d, doy, hr, mi = datetime_converter(t)
    input_file_time_str = 's'+y+doy+hr+mi
    output_file_time_str = y+m+d+'-'+hr+mi
    plot_time_str = y+'-'+m+'-'+d+' '+hr+':'+mi+' UTC'

    #Searching for the data
    dset, data_check = cmi_file_finder(input_loc, input_file_time_str)

    #Checking that we have data to pull
    if data_check==False:
        continue
    
    #Getting the cmi data, with or without the sza correction
    if sza_switch=='sza':
        cmi = sza.sza_calculator_v2_exact(dset)
    elif sza_switch=='raw':
        cmi = dset.variables['CMI'][:]
    else:
        print ('SZA switch error: Default to raw')
        cmi = dset.variables['CMI'][:]
        sza_switch='raw'

    #Getting additional data for plotting
    geo_crs, extent = plot_data(dset)
    orbital_slot = dset.orbital_slot
    scene_id = dset.scene_id
    platform_id = dset.platform_ID
    output_loc_full = output_loc+platform_id+'-'+scene_id+'/'

    if not os.path.exists(output_loc_full):
        os.makedirs(output_loc_full)


    #Making the plot
    fig = plt.figure(figsize=(14,9))
    ax = fig.add_subplot(1, 1, 1, projection=geo_crs)
    ax.imshow(np.sqrt(cmi), extent=extent, transform=geo_crs, cmap='grey',vmin=0, vmax=1.3)

    ax.coastlines(resolution='10m', color='black', linewidth=1)
    ax.add_feature(ccrs.cartopy.feature.STATES)
    ax.set_extent(pbounds)
    ax.set_title(pltname+'\n'+orbital_slot+' '+scene_id+' - '+plot_time_str)
    
    plt.savefig(output_loc_full+platform_id+'-'+scene_id+'-CH'+ch_num+'-'+sza_switch+'-'+output_file_time_str+'.png', facecolor='w', dpi=300)
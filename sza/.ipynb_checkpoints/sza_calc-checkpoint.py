#==============================================
# A series of functions for calculating the solar zenith angle (SZA) and applying it to GOES imagery
#
# Author: Kevin Thiel (kevin.thiel@ou.edu)
# Created: February 2026
#==============================================

#Importing dependencies for calcualting SZA
from pyorbital.astronomy import cos_zen as pyob_cos_zen #Specific calculation needed for solar zenith angle
import numpy as np
import netCDF4 as nc
import pandas as pd
from datetime import datetime



# Source: https://www.star.nesdis.noaa.gov/atmospheric-composition-training/python_abi_lat_lon.php
# Please acknowledge the NOAA/NESDIS/STAR Aerosols and Atmospheric Composition Science Team if using any of this code in your work/research!
def calculate_degrees(file_id):
    
    # Read in GOES ABI fixed grid projection variables and constants
    x_coordinate_1d = file_id.variables['x'][:]  # E/W scanning angle in radians
    y_coordinate_1d = file_id.variables['y'][:]  # N/S elevation angle in radians
    projection_info = file_id.variables['goes_imager_projection']
    lon_origin = projection_info.longitude_of_projection_origin
    H = projection_info.perspective_point_height+projection_info.semi_major_axis
    r_eq = projection_info.semi_major_axis
    r_pol = projection_info.semi_minor_axis
    
    # Create 2D coordinate matrices from 1D coordinate vectors
    x_coordinate_2d, y_coordinate_2d = np.meshgrid(x_coordinate_1d, y_coordinate_1d)
    
    # Equations to calculate latitude and longitude
    lambda_0 = (lon_origin*np.pi)/180.0  
    a_var = np.power(np.sin(x_coordinate_2d),2.0) + (np.power(np.cos(x_coordinate_2d),2.0)*(np.power(np.cos(y_coordinate_2d),2.0)+(((r_eq*r_eq)/(r_pol*r_pol))*np.power(np.sin(y_coordinate_2d),2.0))))
    b_var = -2.0*H*np.cos(x_coordinate_2d)*np.cos(y_coordinate_2d)
    c_var = (H**2.0)-(r_eq**2.0)
    r_s = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)
    s_x = r_s*np.cos(x_coordinate_2d)*np.cos(y_coordinate_2d)
    s_y = - r_s*np.sin(x_coordinate_2d)
    s_z = r_s*np.cos(x_coordinate_2d)*np.sin(y_coordinate_2d)
    
    # Ignore numpy errors for sqrt of negative number; occurs for GOES-16 ABI CONUS sector data
    np.seterr(all='ignore')
    
    abi_lat = (180.0/np.pi)*(np.arctan(((r_eq*r_eq)/(r_pol*r_pol))*((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
    abi_lon = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)
    
    return abi_lat, abi_lon, x_coordinate_2d, y_coordinate_2d


# First attempt to calculate the solar zenith angle for CMI imagery
# This method calculates sza for each pixel individually
def sza_calculator_v1_exact(dset, sza_threshold=88.85):
    '''
    Inputs:
        dset: netCDF4 dataset
        sza_threshold: Maximum solar zenith angle to apply adjustment (degrees). Default=88.85 degrees.
    
    Outputs:
        sza_cmi: CMI adjusted for its solar zenith angle
    '''

    #Calculating the lats/lons from the geostationary projection coordinates (radians)
    tstart = datetime.now()
    lats, lons, x_coordinate_2d, y_coordinate_2d = calculate_degrees(dset)
    print ('Rad2Deg Calculation: '+str(datetime.now()-tstart))

    #Getting the time to use in solar zenith angle calculation
    time = pd.to_datetime(dset.time_coverage_start)

    tstart = datetime.now()
    #Calculating the cosine of solar zenith angle
    cos_zen_grid = pyob_cos_zen(time, lons, lats)
    print ('cos SZA Calculation: '+str(datetime.now()-tstart))

    #Masking out pixels with a sun angle < 0 and too close to sunrise/sunset
    cos_sza_threshold = np.cos(np.radians(sza_threshold))
    cos_zen_grid[cos_zen_grid<cos_sza_threshold] = 1
    
    #Applying the solar zenith angle adjustment to the ABI reflectance factor data
    cmi = dset['CMI'][:] #Reflectance factor
    sza_cmi = cmi / cos_zen_grid

    return sza_cmi
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
from pyproj import Proj
from pyresample import SwathDefinition, kd_tree



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
    
    return abi_lat, abi_lon


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
    lats, lons = calculate_degrees(dset)
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

# Cleaning up code with a way to load npz data files for lat and lon locations (masked array)
def npz_loader(file_name):
    with np.load(file_name) as npz:
        reconstructed_arr = np.ma.masked_array(**npz)  
    return reconstructed_arr

#Loader for latlon data
def lat_lon_loader_conus(dset):
    if (dset.orbital_slot=='GOES-East')&(dset.spatial_resolution=='0.5km at nadir'):
        lats = npz_loader('conus_latlon/east-conus-05-lat.npz')
        lons = npz_loader('conus_latlon/east-conus-05-lon.npz')
    elif (dset.orbital_slot=='GOES-East')&(dset.spatial_resolution=='1km at nadir'):
        lats = npz_loader('conus_latlon/east-conus-10-lat.npz')
        lons = npz_loader('conus_latlon/east-conus-10-lon.npz')
    elif (dset.orbital_slot=='GOES-East')&(dset.spatial_resolution=='2km at nadir'):
        lats = npz_loader('conus_latlon/east-conus-20-lat.npz')
        lons = npz_loader('conus_latlon/east-conus-20-lon.npz')
    elif (dset.orbital_slot=='GOES-West')&(dset.spatial_resolution=='0.5km at nadir'):
        lats = npz_loader('conus_latlon/west-conus-05-lat.npz')
        lons = npz_loader('conus_latlon/west-conus-05-lon.npz')
    elif (dset.orbital_slot=='GOES-West')&(dset.spatial_resolution=='1km at nadir'):
        lats = npz_loader('conus_latlon/west-conus-10-lat.npz')
        lons = npz_loader('conus_latlon/west-conus-10-lon.npz')
    elif (dset.orbital_slot=='GOES-West')&(dset.spatial_resolution=='2km at nadir'):
        lats = npz_loader('conus_latlon/west-conus-20-lat.npz')
        lons = npz_loader('conus_latlon/west-conus-20-lon.npz')
    
    return lats, lons
        

# Second attempt to calculate the solar zenith angle for CMI imagery faster
# This method calculates sza for each pixel individually ONLY for MESO data (CONUS is pulled from pre-saved npz file)
def sza_calculator_v2_exact(dset, sza_threshold=88.85):
    '''
    Inputs:
        dset: netCDF4 dataset
        sza_threshold: Maximum solar zenith angle to apply adjustment (degrees). Default=88.85 degrees.
    
    Outputs:
        sza_cmi: CMI adjusted for its solar zenith angle
    '''

    if dset.scene_id=='CONUS':
        tstart = datetime.now()
        lats, lons = lat_lon_loader_conus(dset)
        print ('Data Load: '+str(datetime.now()-tstart))
    else:
        #Calculating the lats/lons from the geostationary projection coordinates (radians)
        tstart = datetime.now()
        lats, lons = calculate_degrees(dset)
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

#============================================
# RGB functions
#============================================

#Convert over to only positive lons for resampling
def neg_lons(lon):
    for iy, ix in np.ndindex(lon.shape):
        if lon[iy,ix]>180:
            lon[iy,ix] = lon[iy,ix]-360
    return lon

#Resampling code
def resample(input_field, input_lats, input_lons, target_lats, target_lons):

    #Adjusting the negative longitudes on the ERA5 Grid
    #Creating 2D grids of the ERA data
    input_lons = neg_lons(input_lons)

    #Creating the swath definitions
    target_swath = SwathDefinition(lons=target_lons,lats=target_lats)
    input_swath = SwathDefinition(lons=input_lons,lats=input_lats)

    #Resampling using a KD-tree to fit the data to a grid
    output = kd_tree.resample_nearest(source_geo_def=input_swath,
                            data=input_field,
                            target_geo_def=target_swath,
                            radius_of_influence=4e4)
    
    return output

#Quick function to convert simple values into color space
def color_space_converter(val, vmin, vmax, vgamma):
    color_val = ((val-vmin) / (vmax-vmin))**(1/vgamma)
    return color_val

#Function for creating the day cloud phase distinction RGB
def rgb_creator_dcpd(r_dset,g_dset,b_dset, scene, sza_adjustment=False):
    '''
    Creates the Day Cloud Phase Distinction RGB (default) from the three RGB bands as input
    
    :param r_dset: red band dataset (netCDF)
    :param g_dset: green band dataset (netCDF)
    :param b_dset: blue band dataset (netCDF)
    :param scene: 'CONUS' or 'MESO' 

    :returns: rgb
    '''
    #Constants
    rmin = 280.65
    rmax = 219.65
    rgamma = 1.
    gmax = 0.78
    gmin = 0.0
    ggamma = 1.
    bmax = 0.59
    bmin = 0.01
    bgamma = 1.

    #Reading in the datsets to get the values
    r = r_dset.variables['CMI'][:]
    if sza_adjustment==True:
        g = sza_calculator_v2_exact(g_dset)
        b = sza_calculator_v2_exact(b_dset)
    else:
        g = g_dset.variables['CMI'][:]
        b = b_dset.variables['CMI'][:]

    #Getting the lat/lon values for each scene
    if scene=='CONUS':
        r_lats, r_lons = lat_lon_loader_conus(r_dset)
        g_lats, g_lons = lat_lon_loader_conus(g_dset)
        b_lats, b_lons = lat_lon_loader_conus(b_dset)
    else:
        r_lats, r_lons = calculate_degrees(r_dset)
        g_lats, g_lons = calculate_degrees(g_dset)
        b_lats, b_lons = calculate_degrees(b_dset)
    
    #Resampling so all products have the same 0.5 km resolution as the green band
    r = resample(r, r_lats, r_lons, g_lats, g_lons)
    b = resample(b, b_lats, b_lons, g_lats, g_lons)

    #Converting values from the bands into the RGB colorspace
    R = color_space_converter(r, rmin, rmax, rgamma)
    G = color_space_converter(g, gmin, gmax, ggamma)
    B = color_space_converter(b, bmin, bmax, bgamma)

    #Stacking all three together to plot as one RGB
    rgb = np.dstack([R,G,B])

    return rgb
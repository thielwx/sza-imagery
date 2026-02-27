#==============================================
# A script that generates the conus_latlon files at different resolutions (0.5, 1, 2)
#
# Author: Kevin Thiel (kevin.thiel@ou.edu)
# Created: February 2026
#==============================================

import netCDF4 as nc
import numpy as np
import os
import sza_calc as sza

#Location of reference files automatically uploaded to Github and used as a refernece for calculating the conus latlon files
file_loc = '../ref_files'
output_loc = 'conus_latlon/'

if not os.path.exists(output_loc):
    os.makedirs(output_loc)

east_conus_05_fstring = 'OR_ABI-L2-CMIPC-M6C02_G19_s20260580601180_e20260580603553_c20260580604055.nc'
east_conus_10_fstring = 'OR_ABI-L2-CMIPC-M6C05_G19_s20260580601180_e20260580603553_c20260580604021.nc'
east_conus_20_fstring = 'OR_ABI-L2-CMIPC-M6C13_G19_s20260580601180_e20260580603564_c20260580604061.nc'
west_conus_05_fstring = 'OR_ABI-L2-CMIPC-M6C02_G18_s20260580601187_e20260580603561_c20260580604045.nc'
west_conus_10_fstring = 'OR_ABI-L2-CMIPC-M6C05_G18_s20260580601187_e20260580603561_c20260580604065.nc'
west_conus_20_fstring = 'OR_ABI-L2-CMIPC-M6C13_G18_s20260580601187_e20260580603572_c20260580604034.nc'

file_list = [file_loc+east_conus_05_fstring,
             file_loc+east_conus_10_fstring,
             file_loc+east_conus_20_fstring,
             file_loc+west_conus_05_fstring,
             file_loc+west_conus_10_fstring,
             file_loc+west_conus_20_fstring]

for f in file_list:
    #Reading in the file
    dset = nc.Dataset(f, 'r')

    #Calculating the sza
    lat, lon = sza.calculate_degrees(dset)

    #Getting the extra info for saving out the file
    orbital_slot = dset.orbital_slot
    sat_resolution = dset.spatial_resolution

    if orbital_slot == 'GOES-East':
        slot = 'east'
    elif orbital_slot == 'GOES-West':
        slot = 'west'

    if sat_resolution == '0.5km at nadir':
        res = '05'
    elif sat_resolution == '1km at nadir':
        res = '10'
    elif sat_resolution == '2km at nadir':
        res = '20'

    #Saving out the file as compressed npz files
    lat_str = output_loc+slot+'-conus-'+res+'-lat.npz'
    lon_str = output_loc+slot+'-conus-'+res+'-lon.npz'

    np.savez_compressed(lat_str, data=lat.data, mask=lat.mask)
    np.savez_compressed(, data=lon.data, mask=lon.mask)

    print ('Files Saved:')
    print (lat_str)
    print (lon_str)
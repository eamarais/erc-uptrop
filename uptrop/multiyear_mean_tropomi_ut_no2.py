#!/usr/bin/python

"""
Brief python script to calculate multiyear seasonal means of TROPOMI UT NO2.

This is done by weighting single year seasonal means by the number of 
observations.

"""
# Import packages:
import glob
import sys
import os
from os import path
from netCDF4 import Dataset

import numpy as np

# Steps to follow:
# (1) Define input:
# Season of interest (eventually change to argparse)
SEAS = "mam"
# Data product of interest:
PRODUCT = "fresco-wide"
# Define altitude range:
ALT_RANGE = '180-450hPa'
# Data directory:
IODIR = "/home/ucfaea1/python/cloud_slice_pal_no2/Data"
# Years of interest:
YEARS = ['2019', '2020', '2021']
# Year range for output:
YrRange = '2019-2021'

# (2) Get files of interest:
tfile1 = path.join(IODIR, 'tropomi-ut-no2-' + PRODUCT + '-07-1x1-' + SEAS + '-' + YEARS[0] + '-' + ALT_RANGE + '-v2.nc')
#tfile1 = glob.glob(tfile_glob_string)

tfile2 = path.join(IODIR, 'tropomi-ut-no2-' + PRODUCT + '-07-1x1-' + SEAS + '-' + YEARS[1] + '-' + ALT_RANGE + '-v2.nc')
#tfile2 = glob.glob(tfile_glob_string)

tfile3 = path.join(IODIR, 'tropomi-ut-no2-' + PRODUCT + '-07-1x1-' + SEAS + '-' + YEARS[2] + '-' + ALT_RANGE + '-v2.nc')

utno2_files = [tfile1, tfile2, tfile3]

# Get number of files:
nfiles = len(utno2_files)

# Initialize:
first = 0

# Loop over files:
for f in range(nfiles):

    # Track progress:
    print('===> Reading file: ', utno2_files[f])

    # (3) Read data from files:
    fh = Dataset(utno2_files[f], mode='r')

    # latitude:
    tlat = fh.variables['lat'][:]

    # longitude:
    tlon = fh.variables['lon'][:]

    # TROPOMI UT NO2:
    tutno2 = fh.variables['utno2'][:]

    # TROPOMI UT NO2 error:
    tutno2err = fh.variables['utno2err'][:]

    # Cloud top pressre range:
    tcld_top_p_range = fh.variables['cld_top_p_range'][:]

    # Cloud top pressure ceiling:
    tcld_top_p_ceil = fh.variables['cld_top_p_ceil'][:]

    # Number of observations:
    tnobs = fh.variables['nobs'][:]

    # Close file:
    fh.close()

    if first==0:
        # Get dimensions:
        nx = len(tlon)
        ny = len(tlat)

        # Define arrays:
        my_utno2 = np.zeros((nx, ny))
        my_utno2err = np.zeros((nx, ny))
        my_cld_top_p_range = np.zeros((nx, ny))
        my_cld_top_p_ceil = np.zeros((nx, ny))
        my_nobs = np.zeros((nx, ny))
        my_nyrs = np.zeros((nx, ny))
        
        # Redefine:
        first = 1

    # Add individual years:
    indices = np.where(~np.isnan(tutno2))
    my_utno2[indices] += tutno2[indices] * tnobs[indices]
    my_utno2err[indices] += tutno2err[indices] * tnobs[indices]
    my_cld_top_p_range[indices] += tcld_top_p_range[indices] * tnobs[indices]
    my_cld_top_p_ceil[indices] += tcld_top_p_ceil[indices] * tnobs[indices]
    my_nobs[indices] += tnobs[indices]
    my_nyrs[indices] += 1

# Get multiyear mean:
my_utno2 = np.where(my_nobs!=0, np.divide(my_utno2, my_nobs), np.nan)
my_utno2err = np.where(my_nobs!=0, np.divide(my_utno2err, my_nobs), np.nan)
my_cld_top_p_range = np.where(my_nobs!=0, np.divide(my_cld_top_p_range, my_nobs), np.nan)
my_cld_top_p_ceil = np.where(my_nobs!=0, np.divide(my_cld_top_p_ceil, my_nobs), np.nan)

# (5) Save data:
outfile = path.join(IODIR, 'tropomi-ut-no2-' + PRODUCT + '-07-1x1-' + SEAS + '-multiyr-mean-' + YrRange + '-' + ALT_RANGE + '-v2.nc')

ncout = Dataset(outfile, mode='w',format="NETCDF4")

ncout.createDimension('lat', ny)
ncout.createDimension('lon', nx)

# create longitude axis:
lon = ncout.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longiitude'
lon[:] = tlon

# Create latitude axis:
lat = ncout.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lat[:] = tlat

# Save data: UT NO2 VMR (gno2vmr), UT NO2 error (gerr), No. of data points (gcnt)
utno2 = ncout.createVariable('utno2', np.float32, ('lon', 'lat'))
utno2.units = 'pptv'
utno2.long_name = 'NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
utno2[:] = my_utno2

utno2err = ncout.createVariable('utno2err', np.float32, ('lon', 'lat'))
utno2err.units = 'pptv'
utno2err.long_name = 'Standard error of the NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
utno2err[:] = my_utno2err

utcldrange = ncout.createVariable('cld_top_p_range', np.float32, ('lon', 'lat'))
utcldrange.units = 'hPa'
utcldrange.long_name = 'Gridded mean range in cloud top pressures used to cloud-slice TROPOMI NO2'
utcldrange[:] = my_cld_top_p_range

utcldceil = ncout.createVariable('cld_top_p_ceil', np.float32, ('lon', 'lat'))
utcldceil.units = 'hPa'
utcldceil.long_name = 'Gridded mean ceiling of cloud top pressures used to cloud-slice TROPOMI NO2'
utcldceil[:] = my_cld_top_p_ceil

nobs = ncout.createVariable('nobs', np.float32, ('lon', 'lat'))
nobs.units = 'unitless'
nobs.long_name = 'Number of observations in each gridsquare used to obtain cloud-sliced UT NO2 mixing ratios'
nobs[:] = my_nobs

nyrs = ncout.createVariable('nyrs', np.float32, ('lon', 'lat'))
nyrs.units = 'unitless'
nyrs.long_name = 'Number of years in each grid used to obtain multiyear mean'
nyrs[:] = my_nyrs

ncout.close()




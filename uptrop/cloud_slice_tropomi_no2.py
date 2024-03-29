#!/usr/bin/python

"""
Process and apply the cloud-slicing approach to partial columns of NO2 from S5P/TROPOMI for December 2018 to May 2021.

The default is to obtain seasonal means at 1x1 for partial columns above clouds with cloud fraction >=0.7 and within the cloud top pressure range of 450-180 hPa.

Options are available to use cloud information from either the FRESCO-S or ROCINN-CAL cloud product, to obtain seasonal means at 2x2.5 or 4x5, and to use a cloud fraction threshold of 0.8, 0.9, or 1.0.

.. code-block:: bash

    usage:
           [-h] [--trop_dir TROP_DIR] [--out_dir OUT_DIR] [--season SEASON]
           [--start_date START_DATE] [--end_date END_DATE] [--grid_res GRID_RES]
           [--cloud_product CLOUD_PRODUCT] [--cloud_threshold CLOUD_THRESHOLD]
           [--pmin PMIN] [--pmax PMAX]

    optional arguments:
      -h, --help            show this help message and exit
      --trop_dir TROP_DIR   Directory containing tropomi data
      --out_dir OUT_DIR     Directory to contain finished netcdf4
      --season SEASON       Can be jja, son, djf, mam
      --start_date START_DATE
                            Start date of processing window (yyyy-mm-dd)
      --end_date END_DATE   End date of processing window (yyyy-mm-dd)
      --grid_res GRID_RES   Can be 1x1, 2x25, 4x5
      --cloud_product CLOUD_PRODUCT
                            can be fresco or dlr-ocra
      --cloud_threshold CLOUD_THRESHOLD
                            recommended value is 07. Can also test 08, 09, 10
      --pmin PMIN           Lower bound on cloud height. Defaults to 180.
      --pmax PMAX           Upper bound on cloud height. Defaults to 450.


"""

import glob
import argparse
import sys
import os
from os import path
from netCDF4 import Dataset
import datetime as dt
import re

import numpy as np
#from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import matplotlib.pyplot as plt
from dateutil import rrule as rr

# Import hack
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))

from uptrop.height_pressure_converter import alt2pres
from uptrop.cloud_slice_no2 import cldslice, CLOUD_SLICE_ERROR_ENUM
from uptrop.date_file_utils import season_to_date, get_tropomi_file_list, get_ocra_file_list, get_date


class CloudFileDateMismatch(Exception):
    pass

class UnequalColumnException(Exception):
    pass


class GridAggregator:
    """A class for aggregating higher-resolution data into grid squares"""
    def __init__(self, dellat, dellon,pmin,pmax):
        """Creates a grid aggregator across the entire world of resolution dellon, dellat

        :param dellat: vertical size of the aggregation grid in decimal degrees
        :type dellat: float
        :param dellon: Horizontal size of the aggregation grid in decimal degrees
        :type dellon: float

        :returns: A GridAggregator object
        :rtype: GridAggregator"""
        self.postfilt = []
        self.out_lon = np.arange(-180, 180 + dellon, dellon)
        self.out_lat = np.arange(-90, 90 + dellat, dellat)

        self.xdim = len(self.out_lon)
        self.ydim = len(self.out_lat)

        # Seasonal mean output arrays
        self.gno2vmr = np.zeros((self.xdim, self.ydim))  # NO2 VMR in pptv
        self.gcnt = np.zeros((self.xdim, self.ydim))  # No. of data points (orbits)
        self.gerr = np.zeros((self.xdim, self.ydim))  # Weighted error
        self.gwgt = np.zeros((self.xdim, self.ydim))  # Weights
        self.pcld_range = np.zeros((self.xdim, self.ydim))  # Cloud top pressure range used for cloud-slicing
        self.pcld_ceil = np.zeros((self.xdim, self.ydim))  # Ceiling (min pressure) of the cloud top pressure range

        self.file_count = 0
        self.current_max_points = 0

        # Calculate pressure range and centre for Gaussian weighting:
        self.press_sigma = 0.5*(pmax - pmin)
        self.press_mid   = pmin + self.press_sigma

        # Define cloud pressure difference threshold:
        # Use half the range cloud height range so that code is versatile and
        # because this is consistent with Marais et al. (2021) for the upper troposphere.
        self.diff_cldh_thold = self.press_sigma
        #if ((pmin==180) & (pmax==450)): self.diff_cldh_thold=140
        #if ((pmin==180) & (pmax==320)): self.diff_cldh_thold=100
        #if ((pmin==320) & (pmax==450)): self.diff_cldh_thold=100
        #if ((pmin==280) & (pmax==400)): self.diff_cldh_thold=80

        self.loss_count = {
            "too_few_points": 0,
            "low_cloud_height_range": 0,
            "low_cloud_height_std": 0,
            "large_error": 0,
            "sig_diff_from_zero": 0,
            "no2_outlier": 0,
            "non_uni_strat": 0,
        }

        # Members from add_trop_data_to_gridsquare
        self.gno2 = None
        self.gstrat = None
        self.gcldp = None
        self.cntloop = None

        self.cloud_slice_count = 0

    def initalise_grid(self):
        """Zeros the gno2, gstrat, gcldp and cntloop members with empty appendable lists"""
        self.gno2 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.gstrat = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.gcldp = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.cntloop = [[0 for n in range(self.ydim)] for m in range(self.xdim)]

    def grid_trop_data(self, trop_data):
        """Allocates the strat, gc_data and cloud pressure in trop_data into the gno2, gstrat and gcldp grid

        :param trop_data: an instance of TropomiData
        :type trop_data: uptrop.tropomi_ut.gc_data.TropomiData
        """

        for geo_tot_value, trop_lat, trop_lon, strat_no2_val, cloud_pressure \
         in zip(trop_data.geototvcd, trop_data.lats, trop_data.lons, trop_data.stratno2, trop_data.cldpres):
            # Skip over pixels where total column is less than stratospheric
            # column. This addresses positive bias in the cloud-sliced results
            # at low concentrations of UT NO2:
            if geo_tot_value < strat_no2_val:
                continue
            
            # Find nearest gridsquare in output grid:
            p = np.argmin(abs(self.out_lon - trop_lon))
            q = np.argmin(abs(self.out_lat - trop_lat))

            # Convert columns from mol/m2 to molec/cm2:
            tvcdno2 = np.multiply(geo_tot_value, trop_data.no2sfac)
            tstrat = np.multiply(strat_no2_val, trop_data.no2sfac)

            # Get relevant data in each output grid square:
            self.gno2[p][q].append(tvcdno2)
            self.gstrat[p][q].append(tstrat)
            self.gcldp[p][q].append(cloud_pressure)

            # Increment indices:
            self.cntloop[p][q] += 1

        # Save % valid points retained to print out average at end of routine:
        self.postfilt.append(100. * (trop_data.tcnt / trop_data.inicnt))
        self.file_count += 1

    def apply_cloud_slice(self, n_slices=40):
        """Applies the cloud slicing algorithm to each square in the populated grid.

        This method walks over every gridsquare in gno2/gstrat/gcldp.
        After sanity checks and conversion to hPa, that pixel's data is split into
        [len(pixel)/n_slices] subsets, each of n_slices length. The split is alternating,
        so in the cast of three subsets pixel_data[0] goes to subset 0, pixel_data[1] to subset 1,
        pixel_data[2] to subset 2, pixel_data[3] to subset 0, pixel_data[4] to subset 1 and so on.
        Each subset then has apply_slice applied to it.
        If a pixel has less than 100 data points in it, subsetting is skipped and apply_slice is
        run on the pixel's entire dataset.

        :param n_slices: The number of pixels to apply to each subset
        :type n_slices: int
        """
        # Estimate daily mean VMRs from the clustered data:
        for i in range(self.xdim):
            for j in range(self.ydim):

                tcolno2 = self.gno2[i][j]
                strat = self.gstrat[i][j]
                tcld = self.gcldp[i][j]

                # Only loop over grids with relevant data (identified as
                # non-empty lists):
                if not tcolno2:
                    continue

                # Skip if fewer than 10 points:
                if len(tcolno2) < 10:
                    self.loss_count["too_few_points"] += 1
                    continue

                # Convert from Pa to hPa for intput to the cloud-slicing algorithm:
                tcolno2 = np.multiply(tcolno2, 1e4)
                tcld = np.multiply(tcld, 1e-2)
                strat = np.array(strat)

                # Error check that the sizes of the arrays are equal:
                if (len(tcld) != len(tcolno2)):
                    print('Array sizes ne: cloud height and partial column')
                    raise UnequalColumnException

                # Skip scenes with non-uniform stratosphere using the
                # same threshold as is used for GEOS-Chem:
                if (np.std(strat) / np.mean(strat)) > 0.02:
                    self.loss_count["non_uni_strat"] += 1
                    continue

                # Get number of points:
                npnts = len(tcld)
                if npnts > self.current_max_points:
                    self.current_max_points = npnts
                    print(self.current_max_points, flush=True)

                # Use cloud_slice_ut_no2 function to get NO2 mixing
                # ratio from cloud-slicing:
                # Change this slightly to exploit higher spatial resolution of TROPOMI to
                # increase the number of scenes retrieved:
                if ((npnts >= 10) & (npnts < 50)):
                    self.add_slice(i,j,tcld,tcolno2)

                elif (npnts >= 50):
                    # Define number of iterations:
                    stride = round(npnts / n_slices)
                    nloop = list(range(stride))
                    for w in nloop:
                        subset_t_col_no2 = tcolno2[w::stride]
                        subset_t_cld = tcld[w::stride]
                        self.add_slice(i, j, subset_t_cld, subset_t_col_no2)

    def add_slice(self, i, j, t_cld, t_col_no2):
        """Extracts the upper troposphere gc_data, gc_data error and mean cloud pressure for grid square [i,j]

        This method uses the cloud-slicing function :ref:`uptrop.cloud_slice_ut_no2.cldslice`
        Once calculated, the a weighting is derived from cloud pressure.
        The weighted upper tropospheric gc_data and error is added to the rolling total for this season.
        If the cloud slicing fails, then the reason is added to loss_count for the end report.

        :param i: X-index of grid square
        :type i: int
        :param j: Y-index of grid square
        :type j: int
        :param t_cld: A list of cloud pressures
        :type t_cld: list of floats
        :param t_col_no2: A list of gc_data values, of same length as t_cld
        :type t_col_no2: list of floats
        """
        utmrno2, utmrno2err, stage_reached, mean_cld_pres = cldslice(t_col_no2, t_cld,self.diff_cldh_thold)
        # Calculate weights:
        #gaus_wgt = np.exp((-(mean_cld_pres - 315) ** 2) / (2 * 135 ** 2))
        # Skip if approach didn't work (i.e. cloud-sliced UT NO2 is NaN):
        # Drop out after the reason for data loss is added to loss_count.
        if np.isnan(utmrno2) or np.isnan(utmrno2err):
            # Cloud-slicing led to nan, but due to rma regression rather
            # than data filtering (these are rare):
            if (stage_reached == 0):
                print("Cloud-sliced NO2 NAN for pixel i:{} j:{}".format(i, j), flush=True)
                return
            self.loss_count[CLOUD_SLICE_ERROR_ENUM[stage_reached]] += 1
            # print("Cloud-slice exception {} in pixel i:{} j:{}".format(
            #    CLOUD_SLICE_ERROR_ENUM[stage_reached], i, j))
        else:
            # Calculate weights:
            gaus_wgt = np.exp((-(mean_cld_pres - self.press_mid) ** 2) / (2 * self.press_sigma ** 2))
            self.gno2vmr[i, j] += np.multiply(utmrno2, gaus_wgt)
            self.pcld_range[i, j] += np.nanmax(t_cld) - np.nanmin(t_cld)
            self.pcld_ceil[i, j] += np.nanmin(t_cld)
            self.gwgt[i, j] += gaus_wgt
            self.gerr[i, j] += np.multiply(utmrno2err, gaus_wgt)
            self.gcnt[i, j] += 1
            self.cloud_slice_count += 1

    def calc_seasonal_means(self):
        """Calculates the mean no2 mixing ratio using Gaussian weights or counts.
        This is to be applied at the end of processing to get the final data that will be saved and plotted.
        """
        self.mean_gno2vmr = np.divide(self.gno2vmr, self.gwgt, where=self.gcnt != 0)
        self.mean_gerr = np.divide(self.gerr, self.gwgt, where=self.gcnt != 0)
        self.mean_gwgt = np.divide(self.gwgt, self.gcnt, where=self.gcnt != 0)
        self.mean_cld_p_range = np.divide(self.pcld_range, self.gcnt, where=self.gcnt != 0)
        self.mean_cld_p_ceil = np.divide(self.pcld_ceil, self.gcnt, where=self.gcnt != 0)
        self.mean_gno2vmr[self.gcnt == 0] = np.nan
        self.mean_gerr[self.gcnt == 0] = np.nan
        self.mean_gwgt[self.gcnt == 0] = np.nan
        self.mean_cld_p_range[self.gcnt == 0] = np.nan
        self.mean_cld_p_ceil[self.gcnt == 0] = np.nan
        self.gcnt[self.gcnt == 0] = np.nan   # Watch out for this rewriting of gcnt in the future

    def print_report(self):
        """Prints a report of useful data and reasons for data loss."""
        print('Max no. of data points in a gridsquare: ', np.nanmax(self.gcnt), flush=True)
        # Track reasons for data loss:
        print('(1) Too few points: ', self.loss_count["too_few_points"], flush=True)
        print('(2) Low cloud height range: ', self.loss_count["low_cloud_height_range"], flush=True)
        print('(3) Low cloud height std dev: ', self.loss_count["low_cloud_height_std"], flush=True)
        print('(4) Large error: ', self.loss_count["large_error"], flush=True)
        print('(5) Significantly less than zero: ', self.loss_count["sig_diff_from_zero"], flush=True)
        print('(6) Outlier (NO2 > 200 pptv): ', self.loss_count["no2_outlier"], flush=True)
        print('(7) Non-uniform stratosphere: ', self.loss_count["non_uni_strat"], flush=True)
        print('(8) Successful retrievals: ', self.cloud_slice_count, flush=True)
        print('(9) Total possible points: ', (sum(self.loss_count.values()) + self.cloud_slice_count), flush=True)
        print('Mean % points retained: ', np.mean(self.postfilt),flush=True)

    def save_to_netcdf(self, out_file):
        """Saves the seasonal_means to out_file as a netcdf4

        Call after calling calc_seasonal_means at least once.

        :param out_file: Location to save the netcdf4.
        :type out_file: str"""
        ncout = Dataset(out_file, mode='w',format="NETCDF4")

        ncout.createDimension('lat', self.ydim)
        ncout.createDimension('lon', self.xdim)

        # create longitude axis:
        lon = ncout.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longiitude'
        lon[:] = self.out_lon

        # Create latitude axis:
        lat = ncout.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lat[:] = self.out_lat

        # Save data: UT NO2 VMR (gno2vmr), UT NO2 error (gerr), No. of data points (gcnt)
        utno2 = ncout.createVariable('utno2', np.float32, ('lon', 'lat'))
        utno2.units = 'pptv'
        utno2.long_name = 'NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
        utno2[:] = self.mean_gno2vmr

        utno2err = ncout.createVariable('utno2err', np.float32, ('lon', 'lat'))
        utno2err.units = 'pptv'
        utno2err.long_name = 'Standard error of the NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
        utno2err[:] = self.mean_gerr

        utcldrange = ncout.createVariable('cld_top_p_range', np.float32, ('lon', 'lat'))
        utcldrange.units = 'hPa'
        utcldrange.long_name = 'Gridded mean range in cloud top pressures used to cloud-slice TROPOMI NO2'
        utcldrange[:] = self.mean_cld_p_range

        utcldceil = ncout.createVariable('cld_top_p_ceil', np.float32, ('lon', 'lat'))
        utcldceil.units = 'hPa'
        utcldceil.long_name = 'Gridded mean ceiling of cloud top pressures used to cloud-slice TROPOMI NO2'
        utcldceil[:] = self.mean_cld_p_ceil

        nobs = ncout.createVariable('nobs', np.float32, ('lon', 'lat'))
        nobs.units = 'unitless'
        nobs.long_name = 'Number of observations in each gridsquare used to obtain cloud-sliced UT NO2 mixing ratios'
        nobs[:] = self.gcnt

        ncout.close()

    def plot_data(self, out_file):
        """Plots the seasonal_means to screen."""
        # Plot the data:
        X, Y = np.meshgrid(self.out_lon, self.out_lat, indexing='ij')

        # Define plot parameters:
        nplots = 3
        max_val = [80, 80, 20]
        nbound = [21, 21, 21]
        unit = ['[pptv]','[pptv]','unitless']
        plot_title = ['TROPOMI cloud-sliced NO2',
                      'TROPOMI cloud-sliced NO2 error',
                      'Number of cloud-sliced data points']

        fig, ax = plt.subplots(3, 1, figsize=(5,11), subplot_kw=dict(projection=ccrs.PlateCarree()))

        # Plot the subplots:
        for i in range(nplots):

            # Define what data to plot:
            if i==0: plot_vals = np.squeeze(self.mean_gno2vmr)
            if i==1: plot_vals = np.squeeze(self.mean_gerr)
            if i==2: plot_vals = np.squeeze(self.gcnt)

            if i==0: 
                tickval = [0,25,50,75,100]
            else:
                tickval = [0,5,10,15,20]

            ax[i].coastlines(resolution='50m')

            ax[i].set_extent([-179.9, 179.9, -75, 75], crs=ccrs.PlateCarree())

            data_crs = ccrs.PlateCarree()
            bounds = np.linspace(0, max_val[i],nbound[i])

            c = ax[i].pcolormesh(X,Y, plot_vals, transform=data_crs,cmap='jet', 
                                 vmin=0, vmax=max_val[i])

            cb = fig.colorbar(c, ax=ax[i],label=unit[i], orientation='horizontal',
                              shrink=0.5,pad=0.01,boundaries=bounds,ticks=tickval )
            
            cb.ax.tick_params(labelsize=10, direction='in', length=6)

        #plt.savefig(out_file, format='ps')
        plt.show()

class TropomiData:
    """A class for extracting, preprocessing and containing data from a s5p tropomi file."""

    # Initialize _everything_ before the first read.
    def __init__(self, file_path, no2_prod):
        """Reads the tropomi file at file_path and prepares it for processing.

        :param file_path: Path to the netcdf4 file containing the tropomi data
        :type file_path: str
        """

        self.file_name = path.basename(file_path)
        print('Processing: ', self.file_name, flush=True)
        self.date = get_date(self.file_name)
        self.no2_prod = no2_prod

        # Members straight from trop body
        self.no2sfac = None
        self.qasfac = None
        self.fillval = None
        self.tlons = None
        self.tlats = None
        self.tscdno2 = None
        self.stratno2_og = None
        #self.tscdno2err = None  # preserve for future use
        #self.stratno2err = None # preserve for future use
        self.tstratamf = None
        self.qaval = None
        #self.aai = None
        self.sza = None
        self.vza = None

        # Members from geometric column
        self.tamf_geo = None
        self.tgeotropvcd = None

        # Members from bias correction
        self.tstratno2 = None

        # Members from filtering
        self.inicnt = None
        self.tcnt = None

        # Data members from trimming + reshaping
        self.lons = None
        self.lats = None
        self.cldpres = None # Should this be here?
        self.stratno2 = None
        self.amf_geo = None
        self.geototvcd = None

        # Do initialisation
        self.read_trop_file(file_path)

    def read_trop_file(self, file_path):
        """Reads the data at file_path into this object. Called by the constructor.

        :param file_path: Path to the netcdf4 file containing the tropomi data
        :type file_path: str
        """
        fh = Dataset(file_path, mode='r')
        self.fh = fh

        # no2sfac, qasfac, and fillval only need to be read in once, so could
        # just be read in on the first iteration.
        # NO2 conversion factor:
        no2sfac = fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric' \
                                                 '_column'].multiplication_factor_to_convert_to_molecules_percm2

        # QA flag scale factor:
        qasfac = fh.groups['PRODUCT'].variables['qa_value'].scale_factor

        # NO2 fill/missing value:
        fillval = fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric' \
                                                 '_column']._FillValue
        # Extract data of interest:

        # Geolocation data:
        glons = fh.groups['PRODUCT'].variables['longitude'][:]
        tlons = glons[0, :, :]
        glats = fh.groups['PRODUCT'].variables['latitude'][:]
        tlats = glats[0, :, :]

        # Column data:
        # (1) Total slant column:
        gscdno2 = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                      variables['nitrogendioxide_slant_column_density'][:]
        tscdno2 = gscdno2.data[0, :, :]
        # (2) Stratospheric vertical column:
        gstratno2 = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                        variables['nitrogendioxide_stratospheric_column'][:]
        stratno2_og = gstratno2.data[0, :, :]

        # Precisions/Uncertainties of column data (preserve for future use)
        # (1) Total slant column uncertainty:
        #gscdno2err = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
        #                 variables['nitrogendioxide_slant_column_density_precision'][:]
        #tscdno2err = gscdno2err.data[0, :, :]
        # (2) Stratospheric vertical column uncertainty:
        #stratno2err = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
        #                  variables['nitrogendioxide_stratospheric_column_precision'][0, :, :]
        # Stratospheric AMF:
        gstratamf = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                        variables['air_mass_factor_stratosphere'][:]
        tstratamf = gstratamf.data[0, :, :]

        # QA value:
        qaval = fh.groups['PRODUCT'].variables['qa_value'][0, :, :]

        # Aerosol absorbing index:
        # (Preserving for future use)
        #taai = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
        #           variables['aerosol_index_354_388'][:]
        #aai = taai.data[0, :, :]
        #aai = np.where(aai > 1e30, np.nan, aai)

        # Solar zenith angle (degrees):
        tsza = fh.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']. \
                   variables['solar_zenith_angle'][:]
        sza = np.ma.getdata(tsza[0, :, :])

        # Viewing zenith angle (degrees):
        tvza = fh.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']. \
                   variables['viewing_zenith_angle'][:]
        vza = np.ma.getdata(tvza[0, :, :])

        # Setting members
        self.no2sfac = no2sfac
        self.qasfac = qasfac
        self.fillval = fillval
        self.tlons = tlons
        self.tlats = tlats
        self.tscdno2 = tscdno2
        self.stratno2_og = stratno2_og
        #self.tscdno2err = tscdno2err   # Preserve for future use
        #self.stratno2err = stratno2err # Preserve for future use
        self.tstratamf = tstratamf
        self.qaval = qaval
        #self.aai = aai  # Preserve for future use
        self.sza = sza
        self.vza = vza

    def calc_geo_column(self):
        """Calculates the geometric AMF and tropospheric vertical column of this data.
        Includes application of a bias correction to the stratospheric and tropospheric columns. These are obtaind from comparing TROPOMI to Pandora surface observations. 
           The correction addresses an underestimate in TROPOMI stratospheric NO2 variance and a factor of 2 overestimate in TROPOMI tropospheric NO2.
           """

        # Define bias correction values that depend on whether using OFFL or PAL_ NO2 product:
        # No longer needed, as the updated regression approach (Theil Slope) addresses positive bias
        # in background NO2 over remote locations that seemed to necessitate this.
        # It was applied to resolve difference between ground-based (Pandora, MAX-DOAS) and
        # TROPOMI tropospheric columns, but this discrepancy may stem from other issues:
        if self.no2_prod=='PAL':
            trop_correc_fac = 1.0
            
        if self.no2_prod=='OFFL':
            trop_correc_fac = 1.0

        # Calculate the geometric AMF for calculating the tropospheric column:
        tamf_geo = np.add((np.reciprocal(np.cos(np.deg2rad(self.sza)))),
                          (np.reciprocal(np.cos(np.deg2rad(self.vza)))))

        # Bias correct stratospheric column variance:
        if self.no2_prod=='OFFL':
            tstratno2 = np.where(self.stratno2_og != self.fillval, ( (2.5e15 / self.no2sfac) + (self.stratno2_og / 0.87) - (2.8e15 / self.no2sfac)), self.fillval )
        if self.no2_prod=='PAL':
            tstratno2 = np.where(self.stratno2_og != self.fillval, ( (self.stratno2_og / 0.79) - (6.9e14 / self.no2sfac)), self.fillval )

        # Get VCD under cloudy conditions. This is done as the current
        # tropospheric NO2 VCD product includes influence from the prior
        # below clouds:
        # Calculate the stratospheric slant columns using the bias corrected stratospheric column:
        tscdstrat = np.where(self.stratno2_og != self.fillval, (np.multiply(tstratno2, self.tstratamf)), self.fillval )
        # Calculate the tropospheric slant columns:
        ttropscd = np.where(tscdstrat != self.fillval, (np.subtract(self.tscdno2, tscdstrat)), self.fillval )
        # Calculate the tropospheric vertical column using the geometric AMF:
        tgeotropvcd = np.where(ttropscd != self.fillval, (np.divide(ttropscd, tamf_geo)), self.fillval )

        # Decrease tropospheric column by 50% to correct for overestimate
        # compared to Pandora and MAX-DOAS tropospheric columns at Izana:
        tgeotropvcd = np.where(tgeotropvcd != self.fillval, tgeotropvcd / trop_correc_fac, np.nan )
        # The above bias correction has a null effect on the total column, as
        # it just redistributes the relative contribution of the troposphere
        # and the stratosphere.
        # Calculate the correction to the stratospheric column:
        if self.no2_prod=='OFFL':
            tstratno2 = np.where(self.stratno2_og != self.fillval, ( (2.5e15 / self.no2sfac) + (tstratno2 / 0.87) - (2.8e15 / self.no2sfac)), np.nan)
        if self.no2_prod=='PAL':
            tstratno2 = np.where(self.stratno2_og != self.fillval, ( (tstratno2 / 0.87) - (6.9e14 / self.no2sfac)), np.nan)

        # Calculate the total bias corrected column:
        tgeototvcd = np.add(tgeotropvcd, tstratno2)

        # Calculate updated errors after bias correcting.
        # Determine by scaling it by the relative change in stratospheric vertical
        # colum NO2 after applying a bias correction:
        #tstratno2err = np.where(self.stratno2err != self.fillval, np.multiply(self.stratno2err, np.divide(tstratno2, self.stratno2_og)), np.nan)

        # Setting members
        self.tamf_geo = tamf_geo
        self.tgeotropvcd = tgeotropvcd
        self.tstratno2 = tstratno2
        self.tgeototvcd = tgeototvcd

    def cloud_filter_and_preprocess(self, cloud_data, cldthld, pmax, pmin):
        """Filters this tropomi data using the cloud information in cloud_data

        Removes data where
         - There is no cloud data
         - The fraction of cloud is less than the specified cloud threshold
         - Cloud heights are not in the range pmin-pmax
         - Quality value is greater than 0.45

        :param cloud_data: Instance of CloudData
        :type cloud_data: uptrop.tropomi_ut_no2.CloudData
        :param cldthld: The cloud fraction threshold to be used
        :type cldthld: float
        :param pmax: The maximum cloud height to be used in processing
        :type pmax: float
        :param pmin: The minimum cloud height to be used in processing
        :type pmin: float
        """
        # Do date check. This should have already been addressed, but this
        # is an error check in case:
        if self.date != cloud_data.date:
            print('NO2 file: {}, Cloud file: {}'.format(self.date, cloud_data.date), flush=True)
            print('EXITING: Files are not for the same date!', flush=True)
            raise CloudFileDateMismatch

        tgeototvcd = self.tgeototvcd
        # Filter to only include very cloudy scenes at high altitude
        # No. of valid total NO2 column points before apply filtering:
        self.inicnt = np.count_nonzero(~np.isnan(tgeototvcd))

        # Filter out scenes with cloud fractions (and hence cloud pressures) that are nan:
        tgeototvcd = np.where(np.isnan(cloud_data.cldfrac), np.nan, tgeototvcd)

        # Filter out scenes with cloud fraction < cloud threshold:
        tgeototvcd = np.where(cloud_data.cldfrac < cldthld, np.nan, tgeototvcd)

        # Filter out scenes with cloud heights outside the UT range of interest (180-450 hPa):
        tgeototvcd = np.where(cloud_data.tcldpres > pmax * 1e2, np.nan, tgeototvcd)
        tgeototvcd = np.where(cloud_data.tcldpres < pmin * 1e2, np.nan, tgeototvcd)

        # Filter out low quality data (0.45 threshold suggested TROPOMI NO2
        # PI Henk Eskes from the KNMI:
        tgeototvcd = np.where(self.qaval < 0.45, np.nan, tgeototvcd)

        # Filter out scenes with AAI > 1 (could be misclassified as clouds)
        #tgeototvcd = np.where(self.aai > 1., np.nan, tgeototvcd)
        self.tgeototvcd = tgeototvcd

        # No. of points retained after filtering:
        self.tcnt = np.count_nonzero(~np.isnan(self.tgeototvcd))

        # Trim the data to include only those relevant:
        # This also reshapes the data from a 2D to a 1D array:
        self.lons = self.tlons[~np.isnan(tgeototvcd)]
        self.lats = self.tlats[~np.isnan(tgeototvcd)]
        self.stratno2 = self.tstratno2[~np.isnan(tgeototvcd)]
        self.amf_geo = self.tamf_geo[~np.isnan(tgeototvcd)]
        self.geototvcd = self.tgeototvcd[~np.isnan(tgeototvcd)]

        self.cldpres = cloud_data.tcldpres[~np.isnan(tgeototvcd)]

class CloudData:
    """Class for containing the data for cloud filtering and analysis."""
    def __init__(self, file_path, fresco_440, no2_prod, data_type):
        """Reads either the tropomi file (if data_type = 'fresco') or the ocra file (if data_type = 'dlr-ocra')
        at file_path and returns an instance of CloudData. Calls either read_fresco_file or read_ocra_file.

        :param file_path: Path to the file containing cloud_data
        :type file_path: str
        :param data_type: Can be 'fresco' or 'dlr-ocra'
        :type data_type: str
        """
        # Set from file_path
        self.file_name = path.basename(file_path)
        self.date = get_date(self.file_name)  # Assuming for now that ocra and S5P have the same timestamping

        # Set from input:
        self.fresco_440 = fresco_440

        # Set from file contents
        self.cldfrac = None
        self.tcldpres = None
        self.tsnow = None

        # Get cloud fields: 
        if data_type == 'fresco-wide' and no2_prod == 'PAL':
            self.get_fresco_cloud_fields(file_path)
            # Fix snow/ice/coast flag issue in PAL_ files:
            self.fix_snow()
        elif data_type == 'o22cld' and no2_prod == 'PAL':
            self.get_o22cld_cloud_fields(file_path)
            # Fix snow/ice/coast flag issue in PAL_ files:
            self.fix_snow()
        elif data_type == 'dlr-ocra' and no2_prod == 'OFFL':
            # Initialize:
            self.data_parity = True
            self.get_ocra_cloud_fields(file_path)
            self.check_parity()        

    def get_fresco_cloud_fields(self, file_path):
        """Reads and filters the fresco data in a s5p file

        Pixels are dropped if
         - they have type 255 (ocean)
         - they have type 252 (coastline; ATBD treats as 'suspect')
         - they have less thant 1% cover
         - there is a potential misclassification of snow/ice as cloud

        :param file_path: The path to the fresco file
        :type file_path: str
        """
        fh = Dataset(file_path)
        self.fh = fh
        # Cloud fraction (determine whether to use A-band or 440 nm
        #                 effective cloud fraction product):
        #if ( self.fresco_440 ):
        #    tcldfrac = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
        #        variables['cloud_fraction_crb_nitrogendioxide_window'][:]
        #else:
        tcldfrac = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['FRESCO']. \
                   variables['fresco_cloud_fraction_crb'][:]
        self.cldfrac = tcldfrac.data[0, :, :]
        # Cloud top pressure:
        gcldpres = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['FRESCO']. \
                   variables['fresco_cloud_pressure_crb'][:]
        self.tcldpres = gcldpres[0, :, :]
        # Get scene and surface pressure to diagnose clouds misclassified as snow/ice:
        # Apparent scene pressure:
        gscenep = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['FRESCO']. \
                      variables['fresco_apparent_scene_pressure'][:]
        self.tscenep = gscenep.data[0, :, :]

    def get_ocra_cloud_fields(self, file_path):
        """Reads CLOUD_OFFL files with ROCINN-CAL cloud information"""
        # Read data:
        fd = Dataset(file_path)

        # Cloud fraction:
        tcldfrac = fd.groups['PRODUCT'].variables['cloud_fraction'][:]
        cldfrac = tcldfrac.data[0, :, :]

        # Get fill value:
        fillval=(np.max(np.ma.getdata(tcldfrac)))
        if ( fillval<1e35 or fillval==np.nan ):
            print('this method of defining the fill value for dlr-ocra does not work. FIX!!!', flush=True)

        # Cloud top height (m):
        gcldhgt = fd.groups['PRODUCT'].variables['cloud_top_height'][:]
        tcldhgt = np.ma.getdata(gcldhgt[0, :, :])
            
        # Define pressure array of zeros:
        tcldpres = np.zeros(tcldhgt.shape)
            
        # Calculate pressure assuming dry atmosphere using external
        # conversion code (height_pressure_converter.py). There's a cloud
        # top pressure entry in the data file, but this is obtained using
        # ECMWF pressure and might have errors. Diego (DLR cloud product PI
        # recommended I use cloud altitude rather than pressure data):
        hgtind = np.where((tcldhgt != fillval))
        tcldpres[hgtind] = alt2pres(tcldhgt[hgtind])

        # QA value:
        cldqa = fd.groups['PRODUCT'].variables['qa_value'][0, :, :]

        # Snow/ice flag (combined NISE and climatology, so misclassification
        # issues in FRESCO cloud product addressed):
        gsnow = fd.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
            variables['snow_ice_flag'][:]
        self.tsnow = gsnow.data[0, :, :]

        # Set clouds over snow/ice scenes to nan:
        cldfrac = np.where(self.tsnow != 0, np.nan, cldfrac)
        tcldpres = np.where(self.tsnow != 0, np.nan, tcldpres)

        # Set poor quality cloud data to nan:
        self.cldfrac = np.where(cldqa < 0.5, np.nan, cldfrac)
        self.tcldpres = np.where(cldqa < 0.5, np.nan, tcldpres)
        
        # Close file:
        fd.close()

    def get_o22cld_cloud_fields(self, file_path):
        """Reads, filters and preprocesses the data in a dlr-ocra file.

        :ref:`uptrop.height_pressure_converter` is called to convert ocra cloud-top heights to
        pressure values for cross-compatability with fresco data.

        Pixels are dropped if
         - They are over snow/ice scenes
         - They have quality less than 0.5

        :param file_path: The path to the ocra data
        :type file_path: str
        """
        #fh = self.fh
        fh = Dataset(file_path)
        self.fh = fh
        # Cloud fraction:
        tcldfrac = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['O22CLD'].variables['o22cld_cloud_fraction_crb'][:]
        self.cldfrac = tcldfrac.data[0, :, :]
        # Get fill value:
        #fillval=(np.max(np.ma.getdata(tcldfrac)))
        #if ( fillval<1e35 or fillval==np.nan ):
        #    print('this method of defining the fill value for dlr-ocra does not work. FIX!!!', flush=True)
        # Cloud top height (m):
        gcldpres = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['O22CLD'].variables['o22cld_cloud_pressure_crb'][:]
        #tcldhgt = np.ma.getdata(gcldhgt[0, :, :])
        self.tcldpres = gcldpres.data[0, :, :]

        # Apparent scene pressure:
        gscenep = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']['O22CLD']. \
                      variables['o22cld_apparent_scene_pressure'][:]
        self.tscenep = gscenep.data[0, :, :]

       # Define pressure array of zeros:
        #tcldpres = np.zeros(tcldhgt.shape)

        # Calculate pressure assuming dry atmosphere using external
        # conversion code (height_pressure_converter.py). There's a cloud
        # top pressure entry in the data file, but this is obtained using
        # ECMWF pressure and might have errors. Diego (DLR cloud product PI
        # recommended I use cloud altitude rather than pressure data):
        #hgtind = np.where((tcldhgt != fillval))
        #tcldpres[hgtind] = alt2pres(tcldhgt[hgtind])

        # QA value:
        #cldqa = fd.groups['PRODUCT'].variables['qa_value'][0, :, :]

        # Snow/ice flag (combined NISE and climatology, so misclassification
        # issues in FRESCO cloud product addressed):
        #gsnow = fd.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
        #            variables['snow_ice_flag'][:]
        #self.tsnow = gsnow.data[0, :, :]

        # Set clouds over snow/ice scenes to nan:
        #cldfrac = np.where(self.tsnow != 0, np.nan, cldfrac)
        #tcldpres = np.where(self.tsnow != 0, np.nan, tcldpres)

        # Set poor quality cloud data to nan:
        #self.cldfrac = np.where(cldqa < 0.5, np.nan, cldfrac)
        #self.tcldpres = np.where(cldqa < 0.5, np.nan, tcldpres)

        # Close DLR CLOUD file:
        #fd.close()

    def check_parity(self):
        # Skip files if the number of indices are not equal:            
        if self.cldfrac.shape != trop_data.sza.shape:
            print('Cloud product and NO2 indices ne!', flush=True)
            print(self.cldfrac.shape, trop_data.sza.shape, flush=True)
            print('Skipping this swath', flush=True)
            self.data_parity=False
            pass

    def fix_snow(self):
        """Address issue with identifying ocean/coasts/ice/snow

        Pixels dropped if:
         - type 255 (ocean)
         - type 252 (coastline; ATBD treats as 'suspect')
         - < 1% ...  cover
         - potential misclassification of snow/ice as cloud

        """
        fh = self.fh
        # Surface pressure:
        gsurfp = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                     variables['surface_pressure'][:]
        tsurfp = gsurfp.data[0, :, :]

        # Snow/ice flag:
        gsnow = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                    variables['snow_ice_flag'][:]
        tsnow = gsnow.data[0, :, :]
        # Convert all valid snow/ice free flag values (0,255) to 0.
        # Ocean:
        tsnow = np.where(tsnow == 255, 0, tsnow)
        # Coastlines (listed as potential "suspect" in the ATBD document p. 67):
        tsnow = np.where(tsnow == 252, 0, tsnow)
        # Less then 1% snow/ice cover:
        tsnow = np.where(tsnow < 1, 0, tsnow)
        # Snow/ice misclassified as clouds:
        self.tsnow = np.where(((tsnow > 80) & (tsnow < 104) & (self.tscenep > (0.98 * tsurfp))),
                         0, tsnow)
        # Set clouds over snow/ice scenes to nan:
        self.cldfrac = np.where(self.tsnow != 0, np.nan, self.cldfrac)
        self.tcldpres = np.where(self.tsnow != 0, np.nan, self.tcldpres)

        # Close S5P TROPOMI NO2 file:
        fh.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Produces a netCDF and preview plot for upper troposphere NO2")
    parser.add_argument("--trop_dir", help="Directory containing tropomi data")
    parser.add_argument("--out_dir", help="Directory to contain finished netcdf4")
    parser.add_argument("--season", help="Can be jja, son, djf, mam")
    parser.add_argument("--year", help="Can be 2018, 2019, 2020, 2021")
    parser.add_argument("--start_date", help="Start date of processing window (yyyy-mm-dd)")
    parser.add_argument("--end_date", help="End date of processing window (yyyy-mm-dd)")
    parser.add_argument("--grid_res", default='1x1', help="Can be 1x1, 2x25, 4x5")
    parser.add_argument("--cloud_product", default = "fresco-wide", help="can be fresco-wide or o22cld")
    parser.add_argument("--cloud_threshold", default = "07", help="recommended value is 07. Can also test 08, 09, 10")
    parser.add_argument("--fresco_440", type=bool, default = False, help="Use FRESCO-S 440 nm cloud fraction")
    parser.add_argument("--pmin", default=180, type=int, help="Lower bound on cloud height. Defaults to 180.")
    parser.add_argument("--pmax", default=450, type=int, help="Upper bound on cloud height. Defaults to 450.")
    parser.add_argument("--no2_prod", default = "PAL", help="TROPOMI NO2 product name. Can be OFFL or PAL")
    parser.add_argument("--version", default = "v1", help="Version number to append to filename")
    args = parser.parse_args()

    if args.season:
        start_date, end_date = season_to_date(args.season,np.int(args.year))
    else:
        if args.start_date is not None and args.end_date is not None:
            start_date = dt.datetime.strptime(args.start_date, "%Y-%m-%d")
            end_date = dt.datetime.strptime(args.end_date, "%Y-%m-%d")
        else:
            print("Please provide either --season or --start_date and --end_date")
            sys.exit(1)

    if start_date.year == end_date.year:
        yrrange = str(start_date.year)
    else:
        yrrange = str(start_date.year) + "-" + str(end_date.year)

    if args.grid_res == '1x1':
        dellat, dellon = 1, 1
    elif args.grid_res == '2x25':
        dellat, dellon = 2, 2.5
    elif args.grid_res == '4x5':
        dellat, dellon = 4, 5
    elif args.grid_res == '05x05':
        dellat, dellon = 0.5, 0.5
    else:
        print("Invalid grid; values can be 05x05, 1x1, 2x25, 4x5")
        sys.exit(1)

    # Parsing cloud threshold
    cloud_threshold = float(args.cloud_threshold)/10

    # Parsing bollean to use or not use FRESCO-S cloud fraction at 440 nm:
    fresco_440 = args.fresco_440

    date_range = rr.rrule(rr.DAILY, dtstart=start_date, until=end_date)
    trop_files = get_tropomi_file_list(args.trop_dir, date_range, args.no2_prod)
    print('Found total of {} files: '.format(len(trop_files)))
    if args.cloud_product == "fresco-wide":
        cloud_files = trop_files
    elif args.cloud_product == "dlr-ocra":
        cloud_files = get_ocra_file_list(args.trop_dir, date_range)
    elif args.cloud_product == "o22cld":
        cloud_files = trop_files
    else:
        print("Invalid cloud product; can be fresco-wide or o22cld")
        sys.exit(1)

    grid_aggregator = GridAggregator(dellat, dellon,args.pmin,args.pmax)

    # Check for instances where there are cloud files for a TROPOMI swath, but
    # no NO2 files and vice versa:
    # Only do this when using CLOUD_OFFL product:
    if args.cloud_product == "dlr-ocra":
        if len(cloud_files)!=len(trop_files):
            for i, (trop_file, cloud_file) in enumerate(zip(trop_files, cloud_files)):
                # Extract dates for these files:
                trop_date=get_date(trop_file)
                cloud_date=get_date(cloud_file)
                # Check that dates are consistent 
                if ( cloud_date != trop_date ):
                #    if ( cloud_date > trop_date ):
                    # No instances of this found, but preserve in case. Not tested, so commented out:
                    #    # Remove no2 file with no matching cloud file for this
                    #    # swath/orbit:
                    #    del [no2_files[i]]
                    #    # Restart iteration:
                    #    i=0
                    if ( cloud_date < trop_date ):
                        # Remove cloud file with no matching no2 file for this
                        # swath/orbit:
                        del [cloud_files[i]]
                        # Restart iteration:
                        i=0
                        
    for trop_file, cloud_file in zip(trop_files, cloud_files):
        trop_data = TropomiData(trop_file, args.no2_prod)
        cloud_data = CloudData(cloud_file, fresco_440, args.no2_prod, data_type=args.cloud_product)
        if args.cloud_product == "dlr-ocra" and cloud_data.data_parity==False: continue
        trop_data.calc_geo_column()
        trop_data.cloud_filter_and_preprocess(cloud_data, cloud_threshold, args.pmax, args.pmin)
        grid_aggregator.initalise_grid()
        grid_aggregator.grid_trop_data(trop_data)
        grid_aggregator.apply_cloud_slice()
    grid_aggregator.calc_seasonal_means()

    # Define cloud product string for file name:
    str_cld_prod = args.cloud_product
    if ( fresco_440 ):
        str_cld_prod = 'fresco-440nm'
    #out_plot_file = ' '
    #grid_aggregator.plot_data(out_plot_file)

    # Define output file names:
    # Data file and plot file names if season is or isn't specified in the input arguments:
    if args.season is not None: 
        out_data_file = 'tropomi-ut-no2-' + str_cld_prod + '-' + args.cloud_threshold + '-' + args.grid_res + '-' + args.season + '-' + yrrange + '-' + str(args.pmin) + '-' + str(args.pmax) + 'hPa' + '-' + args.version + '.nc'
        out_plot_file = 'tropomi-ut-no2-' + str_cld_prod + '-' + args.cloud_threshold + '-' + args.grid_res + '-' + args.season[0:3] + '-' + yrrange + '-' + str(args.pmin) + '-' + str(args.pmax) + 'hPa' + '-' + args.version + '.ps'
    else:
        out_data_file = 'tropomi-ut-no2-' + str_cld_prod + '-' + args.cloud_threshold + '-' + args.grid_res + '-' + yrrange + '-' + str(args.pmin) + '-' + str(args.pmax) + 'hPa' + '-' + args.version + '.nc'
        out_plot_file = 'tropomi-ut-no2-' + str_cld_prod + '-' + args.cloud_threshold + '-' + args.grid_res + '-' + yrrange + '-' + str(args.pmin) + '-' + str(args.pmax) + 'hPa' + '-' + args.version + '.ps'
    # Complete data file path:
    out_data_file_path = path.join(args.out_dir, 'Data/' + out_data_file)
    #out_data_file = glob.glob(out_data_file_path)
    print("Saving data to: {}".format(out_data_file), flush=True)
    
    # (2) Plot file:
    #out_plot_file_path = path.join('/Images/', out_plot_file)
    out_plot_file_path = path.join(args.out_dir, 'Images/', out_plot_file)
    print("Saving image to: {}".format(out_plot_file_path), flush=True)

    grid_aggregator.print_report()
    grid_aggregator.save_to_netcdf(out_data_file_path)
    #grid_aggregator.plot_data(out_plot_file_path)



    

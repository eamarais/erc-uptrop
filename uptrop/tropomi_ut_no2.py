"""
Process and apply the cloud-slicing approach to partial columns of NO2 from S5P/TROPOMI for June 2019 to May 2020.

The default is to obtain seasonal means at 1x1 for partial columns above clouds with cloud fraction >=0.7 and within the cloud top pressure range of 450-180 hPa.

Options are available to use cloud information from either the FRESCO-S or ROCINN-CAL cloud product, to obtain seasonal means at 2x2.5 or 4x5, and to use a cloud fraction threshold of 0.8, 0.9, or 1.0.

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
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from dateutil import rrule as rr

# Import hack
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..'))

from uptrop.convert_height_to_press import alt2pres
from uptrop.cloud_slice_ut_no2 import cldslice, CLOUD_SLICE_ERROR_ENUM


class CloudFileDateMismatch(Exception):
    pass


class CloudFileShapeMismatch(Exception):
    """
    Raised when FRESCO and DLR files are not the same shape
    """


class UnequalColumnException(Exception):
    pass


class GridAggregator:
    """A class for aggregating higher-resolution data into grid squares"""
    def __init__(self, dellat, dellon):
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

        self.file_count = 0
        self.current_max_points = 0

        self.loss_count = {
            "too_few_points": 0,
            "low_cloud_height_range": 0,
            "low_cloud_height_std": 0,
            "large_error": 0,
            "much_less_than_zero": 0,
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
        """Allocates the strat, no2 and cloud pressure in trop_data into the gno2, gstrat and gcldp grid

        :param trop_data: an instance of TropomiData
        :type trop_data: uptrop.tropomi_ut.no2.TropomiData
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
                if ((npnts >= 10) & (npnts < 100)):
                    self.add_slice(i,j,tcld,tcolno2)

                elif (npnts >= 100):
                    # Define number of iterations:
                    stride = round(npnts / n_slices)
                    nloop = list(range(stride))
                    for w in nloop:
                        subset_t_col_no2 = tcolno2[w::stride]
                        subset_t_cld = tcld[w::stride]
                        self.add_slice(i, j, subset_t_cld, subset_t_col_no2)

    def add_slice(self, i, j, t_cld, t_col_no2):
        """Extracts the upper troposphere no2, no2 error and mean cloud pressure for grid square [i,j]

        This method uses the cloud-slicing function [ref]cld
        Once calculated, the a weighting is derived from cloud pressure.
        The weighted upper tropospheric no2 and error is added to the rolling total for this season.
        If the cloud slicing fails, then the reason is added to loss_count for the end report.

        :param i: X-index of grid square
        :type i: int
        :param j: Y-index of grid square
        :type j: int
        :param t_cld: A list of cloud pressures
        :type t_cld: list of floats
        :param t_col_no2: A list of no2 values, of same length as t_cld
        :type t_col_no2: list of floats
        """
        utmrno2, utmrno2err, stage_reached, mean_cld_pres = cldslice(t_col_no2, t_cld)
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
            gaus_wgt = np.exp((-(mean_cld_pres - 315) ** 2) / (2 * 135 ** 2))
            self.gno2vmr[i, j] += np.multiply(utmrno2, gaus_wgt)
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
        self.mean_gno2vmr[self.gcnt == 0] = np.nan
        self.mean_gerr[self.gcnt == 0] = np.nan
        self.mean_gwgt[self.gcnt == 0] = np.nan
        self.gcnt[self.gcnt == 0] = np.nan   # Watch out for this rewriting of gcnt in the future

    def print_report(self):
        """Prints a report of useful data and reasons for data loss."""
        print('Max no. of data points in a gridsquare: ', np.nanmax(self.gcnt), flush=True)
        # Track reasons for data loss:
        print('(1) Too few points: ', self.loss_count["too_few_points"], flush=True)
        print('(2) Low cloud height range: ', self.loss_count["low_cloud_height_range"], flush=True)
        print('(3) Low cloud height std dev: ', self.loss_count["low_cloud_height_std"], flush=True)
        print('(4) Large error: ', self.loss_count["large_error"], flush=True)
        print('(5) Significantly less then zero: ', self.loss_count["much_less_than_zero"], flush=True)
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

        nobs = ncout.createVariable('nobs', np.float32, ('lon', 'lat'))
        nobs.units = 'unitless'
        nobs.long_name = 'Number of observations in each gridsquare used to obtain cloud-sliced UT NO2 mixing ratios'
        nobs[:] = self.gcnt

        ncout.close()

    def plot_data(self):
        """Plots the seasonal_means to screen."""
        # Plot the data:
        m = Basemap(resolution='l', projection='merc',
                    lat_0=0, lon_0=0, llcrnrlon=-180,
                    llcrnrlat=-75, urcrnrlon=180, urcrnrlat=80)
        X, Y = np.meshgrid(self.out_lon, self.out_lat, indexing='ij')
        xi, yi = m(X, Y)
        plt.subplot(1, 3, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.mean_gno2vmr), vmin=0, vmax=80, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('NO2 VMRs')

        plt.subplot(1, 3, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.mean_gerr), vmin=0, vmax=30, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('NO2 error')

        plt.subplot(1, 3, 3)
        cs = m.pcolor(xi, yi, np.squeeze(self.gcnt), vmin=0., vmax=30, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('Number of points')

        plt.show()


class TropomiData:
    """A class for extracting, preprocessing and containing data from a s5p tropomi file."""

    #This time, I'm initialising _everything_ first before the first read.
    def __init__(self, file_path):
        """Reads the tropomi file at file_path and prepares it for processing.

        :param file_path: Path to the netcdf4 file containing the tropomi data
        :type file_path: str
        """

        self.file_name = path.basename(file_path)
        print('Processing: ', self.file_name, flush=True)
        self.date = get_date(self.file_name)

        # Members straight from trop body
        self.no2sfac = None
        self.qasfac = None
        self.fillval = None
        self.tlons = None
        self.tlats = None
        self.tscdno2 = None
        self.stratno2_og = None
        self.tscdno2err = None
        self.stratno2err = None
        self.tstratamf = None
        self.qaval = None
        self.aai = None
        self.sza = None
        self.vza = None

        # Members from geometric column
        self.tamf_geo = None
        self.tgeotropvcd = None

        # Members from bias correction
        self.tstratno2 = None
        self.tgeototvcd = None
        self.ttropvcd_geo_err = None  # This one doesn't seem to be used

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

        # Precisions/Uncertainties of column data:
        # (1) Total slant column uncertainty:
        gscdno2err = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                         variables['nitrogendioxide_slant_column_density_precision'][:]
        tscdno2err = gscdno2err.data[0, :, :]
        # (2) Stratospheric vertical column uncertainty:
        stratno2err = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                          variables['nitrogendioxide_stratospheric_column_precision'][0, :, :]
        # Stratospheric AMF:
        gstratamf = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                        variables['air_mass_factor_stratosphere'][:]
        tstratamf = gstratamf.data[0, :, :]

        # QA value:
        qaval = fh.groups['PRODUCT'].variables['qa_value'][0, :, :]

        # Aerosol absorbing index:
        taai = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                   variables['aerosol_index_354_388'][:]
        aai = taai.data[0, :, :]
        aai = np.where(aai > 1e30, np.nan, aai)

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
        self.tscdno2err = tscdno2err
        self.stratno2err = stratno2err
        self.tstratamf = tstratamf
        self.qaval = qaval
        self.aai = aai
        self.sza = sza
        self.vza = vza

    def calc_geo_column(self):
        """Calculates the geometric AMF and tropospheric vertical column of this data"""
        # Calculate the geometric AMF:
        tamf_geo = np.add((np.reciprocal(np.cos(np.deg2rad(self.sza)))),
                          (np.reciprocal(np.cos(np.deg2rad(self.vza)))))

        # Get VCD under cloud conditions. This is done as the current
        # tropospheric NO2 VCD product includes influence from the prior
        # below clouds:
        # Calculate the stratospheric slant columns:
        tscdstrat = np.multiply(self.stratno2_og, self.tstratamf)
        # Calculate the tropospheric slant columns:
        ttropscd = np.subtract(self.tscdno2, tscdstrat)
        # Calculate the tropospheric vertical column using the geometric AMF:
        tgeotropvcd = np.divide(ttropscd, tamf_geo)

        # Setting members
        self.tamf_geo = tamf_geo
        self.tgeotropvcd = tgeotropvcd

    def apply_bias_correction(self):
        """Applies bias corrections to this data. These are obtaind from comparing TROPOMI to Pandora surface observations. 
           The correction addresses an underestimate in TROPOMI stratospheric NO2 variance and a factor of 2 overestimate in TROPOMI tropospheric NO2.
           ["""
        # Bias correct stratosphere based on comparison of TROPOMI to Pandora Mauna Loa:
        tstratno2 = np.where(self.stratno2_og != self.fillval,
                             ((self.stratno2_og - (6.6e14 / self.no2sfac)) / 0.86), np.nan)

        # Bias correct troposphere based on comparison of TROPOMI to Pandora Izana:
        tgeotropvcd = np.where(self.tgeotropvcd != self.fillval,
                               self.tgeotropvcd / 2., np.nan)

        # Get the total column as the sum of the bias-corrected components:
        tgeototvcd = np.add(tgeotropvcd, tstratno2)

        # Calculate updated stratospheric NO2 error after bias correcting.
        # Determine by scaling it by the relative change in stratospheric vertical
        # colum NO2 after applying a bias correction:
        tstratno2err = np.where(self.stratno2err != self.fillval,
                                np.multiply(self.stratno2err, np.divide(tstratno2, self.stratno2_og)),
                                np.nan)

        # Calculate error by adding in quadrature individual
        # contributions:
        ttotvcd_geo_err = np.sqrt(np.add(np.square(tstratno2err),
                                         np.square(self.tscdno2err)))
        # Estimate the tropospheric NO2 error as the total error
        # weighted by the relative contribution of the troposphere
        # to the total column, as components that contribute to the
        # error are the same:
        ttropvcd_geo_err = np.multiply(ttotvcd_geo_err,
                                       (np.divide(tgeotropvcd, tgeototvcd)))

        self.tstratno2 = tstratno2
        self.tgeototvcd = tgeototvcd
        self.tgeotropvcd = tgeotropvcd  # Filter applied to member defined in geo_column
        self.ttropvcd_geo_err = ttropvcd_geo_err

    def cloud_filter_and_preprocess(self, cloud_data, cldthld, pmin, pmax):
        """Filters this tropomi data using the cloud information in cloud_data

        Removes data where
         - There is no cloud data
         - The fraction of cloud is less than the specified cloud threshold
         - Cloud heights are not in the range pmin-pmax
         - Quality value is greater than 0.45
         - Aerosol Absorbing Index (AAI) is > 1

        :param cloud_data: Instance of CloudData
        :type cloud_data: uptrop.tropomi_ut_no2.CloudData
        :param cldthld: The cloud fraction threshold to be used
        :type cldthld: float
        :param pmax: The maximum cloud height to be used in processing
        :type pmax: float
        :param pmin: The minimum cloud height to be used in processing
        :type pmin: float
        """
        # Do date check
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
        tgeototvcd = np.where(self.aai > 1., np.nan, tgeototvcd)
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
    def __init__(self, file_path, data_type):
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

        # Set from file contents
        self.cldfrac = None
        self.tcldpres = None
        self.tsnow = None

        # Initialize:
        self.data_parity = True

        if data_type == 'fresco':
            self.read_fresco_file(file_path)
        elif data_type == 'dlr-ocra':
            self.read_ocra_file(file_path)
            self.check_parity()

    def read_fresco_file(self, file_path):
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
        # Cloud fraction:
        tcldfrac = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                       variables['cloud_fraction_crb'][:]
        self.cldfrac = tcldfrac.data[0, :, :]
        # Cloud top pressure:
        gcldpres = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                       variables['cloud_pressure_crb'][:]
        self.tcldpres = gcldpres[0, :, :]

        # Get scene and surface pressure to diagnose clouds misclassified as snow/ice:
        # Apparent scene pressure:
        gscenep = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                      variables['apparent_scene_pressure'][:]
        tscenep = gscenep.data[0, :, :]
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
        self.tsnow = np.where(((tsnow > 80) & (tsnow < 104) & (tscenep > (0.98 * tsurfp))),
                         0, tsnow)
        # Set clouds over snow/ice scenes to nan:
        self.cldfrac = np.where(self.tsnow != 0, np.nan, self.cldfrac)
        self.tcldpres = np.where(self.tsnow != 0, np.nan, self.tcldpres)
        fh.close()

    def read_ocra_file(self, file_path):
        """Reads, filters and preprocesses the data in a dlr-ocra file.

        :ref:`uptrop.convert_height_to_pressure` is called to convert ocra cloud-top heights to
        pressure values for cross-compatability with fresco data.

        Pixels are dropped if
         - They are over snow/ice scenes
         - They have quality less than 0.5

        :param file_path: The path to the ocra data
        :type file_path: str
        """

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
        # conversion code (convert_height_to_press.py). There's a cloud
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

        # Close DLR CLOUD file:
        fd.close()

    def check_parity(self):
        # Skip files if the number of indices are not equal:            
        if self.cldfrac.shape != trop_data.sza.shape:
            print('Cloud product and NO2 indices ne!', flush=True)
            print(self.cldfrac.shape, trop_data.sza.shape, flush=True)
            print('Skipping this swath', flush=True)
            self.data_parity=False

#   TODO: Move these into a seperate file for reuse maybe
def get_tropomi_file_list(trop_dir, date_range):
    """Returns an alphabetically sorted list of Tropomi files
    within a range of dates.

    :param trop_dir: The directory containing the tropomi files
    :type trop_dir: str
    :param date_range: A list of dates. Generation using DateUtil's rrule function is recommended.
    :type date_range: list(datetime)

    :returns: A list of filepaths to tropomi data
    :rtype: list of str
    """
    out = []
    for date in date_range:
        out += (get_tropomi_files_on_day(trop_dir, date))
    return sorted(out)


def get_ocra_file_list(ocra_dir, date_range):
    """Returns an alphabetically sorted list of Ocra files
    within a range of dates.

    :param ocra_dir: The directory containing the ocra files
    :type ocra_dir: str
    :param date_range: A list of dates. Generation using DateUtil's rrule function is recommended.
    :type date_range: list(datetime)

    :returns: A list of filepaths to ocra data
    :rtype: list of str
    """
    out = []
    for date in date_range:
        out += (get_ocra_files_on_day(ocra_dir, date))
    return sorted(out)


def get_tropomi_files_on_day(tomidir, date):
    """Returns a list of tropomi files on a given date.

    Uses the :ref:`get_date` function to extract each candidate file's date from it's filename

    :param tomidir: The directory containing the tropomi files
    :type tomidir: str
    :param date: The date to search for
    :type date: DateTime

    :returns: A list of filepaths to tropomi data
    :rtype: list of str
    """
    # Converts the python date object to a set string representation of time
    # In this case, zero-padded year, month and a datestamp of the Sentinel format
    # See https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
    year = date.strftime(r"%Y")
    month = date.strftime(r"%m")
    datestamp = date.strftime(r"%Y%m%dT")
    tomi_glob_string = path.join(tomidir, 'NO2_OFFL', year, month,'S5P_OFFL_L2__NO2____'+ datestamp + '*')
    tomi_files_on_day = glob.glob(tomi_glob_string)
    print('Found {} tropomi files for {}: '.format(len(tomi_files_on_day), date))
    tomi_files_on_day = sorted(tomi_files_on_day)
    return tomi_files_on_day


def get_ocra_files_on_day(tomidir,date):
    """Returns a list of ocra files on a given date.

    Uses the :ref:`uptrop.tropomi_ut_no2.get_date` function to extract each candidate file's date from it's filename

    :param tomidir: The directory containing the ocra files (usually packaged with tropomi data)
    :type tomidir: str
    :param date: The date to search for
    :type date: DateTime

    :returns: A list of filepaths to ocra data
    :rtype: list of str
    """
    # Get string of day:
    year = date.strftime(r"%Y")
    month = date.strftime(r"%m")
    datestamp = date.strftime(r"%Y%m%dT")
    cld_glob_string = path.join(tomidir, "CLOUD_OFFL", year, month,
                                   'S5P_OFFL_L2__CLOUD__' + datestamp + '*')
    cldfile = glob.glob(cld_glob_string) 
    # Order the files:
    cldfile = sorted(cldfile)
    return cldfile


def get_date(file_name, time_stamp_index = 0):
    """Extracts a datetime object from a filename with a Sentinel timestamp

    See https://regex101.com/r/QNG11l/1 for examples

    :param file_name: The filename to extract the date from
    :type file_name: str
    :param time_stamp_index: Which time-stamp to get the date from if more than one. Defaults to 0.
    :type time_stamp_index: int

    :returns: A DateTime object of the date of the file
    :rtype: DateTime
    """
    # A regular expression that gets Sentinel datestamps out of filenames
    # See https://regex101.com/r/QNG11l/1
    date_regex = r"\d{8}T\d{6}"
    date_string = re.findall(date_regex, file_name)[time_stamp_index]
    # A line for converting Sentinel string reps to datetime
    return dt.datetime.strptime(date_string, r"%Y%m%dT%H%M%S")


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Produces a netCDF and preview plot for upper troposphere NO2")
    parser.add_argument("trop_dir")
    parser.add_argument("out_dir")
    parser.add_argument("--season", default='jja', help="Can be jja, son, djf, mam")
    parser.add_argument("--grid_res", default='1x1', help="Can be 1x1, 2x25, 4x5")
    parser.add_argument("--cloud_product", default = "fresco", help="can be fresco or dlr-ocra")
    parser.add_argument("--cloud_threshold", default = "07", help="recommended value is 07. Can also test 08, 09, 10")
    parser.add_argument("--pmin", default=180, type=int)
    parser.add_argument("--pmax", default=450, type=int)
    args = parser.parse_args()

    if args.season == "jja":
        start_date = dt.datetime(year=2019, month=6, day=1)
        end_date = dt.datetime(year=2019, month=8, day=31)
        yrrange = '2019'
    elif args.season == "son":
        start_date = dt.datetime(year=2019, month=9, day=1)
        end_date = dt.datetime(year=2019, month=11, day=30)
        yrrange = '2019'
    elif args.season == "djf":
        start_date = dt.datetime(year=2019, month=12, day=1)
        end_date = dt.datetime(year=2020, month=2, day=29)  # Beware the leap year here
        yrrange = '2019-2020'
    elif args.season == "mam":
        start_date = dt.datetime(year=2020, month=3, day=1)
        end_date = dt.datetime(year=2020, month=6, day=29)
        yrrange = '2020'
    elif args.season == "test":
        start_date = dt.datetime(year=2020, month=3, day=1)
        end_date = dt.datetime(year=2020, month=3, day=3)
        yrrange = 'TEST'
    else:
        print("Invalid season; can be jja, son, djf, mam")
        sys.exit(1)

    if args.grid_res == '1x1':
        dellat, dellon = 1, 1
    elif args.grid_res == '2x25':
        dellat, dellon = 2, 2.5
    elif args.grid_res == '4x5':
        dellat, dellon = 4, 5
    else:
        print("Invalid grid; values can be 1x1, 2x25, 4x5")
        sys.exit(1)

    # Parsing cloud threshold
    cloud_threshold = float(args.cloud_threshold)/10

    date_range = rr.rrule(rr.DAILY, dtstart=start_date, until=end_date)

    trop_files = get_tropomi_file_list(args.trop_dir, date_range)
    print('Found total of {} files: '.format(len(trop_files)))
    if args.cloud_product == "fresco":
        cloud_files = trop_files
    elif args.cloud_product == "dlr-ocra":
        cloud_files = get_ocra_file_list(args.trop_dir, date_range)
    else:
        print("Invalid cloud product; can be fresco or dlr-ocra")
        sys.exit(1)

    grid_aggregator = GridAggregator(dellat, dellon)

    for trop_file, cloud_file in zip(trop_files, cloud_files):
        trop_data = TropomiData(trop_file)
        cloud_data = CloudData(cloud_file, data_type=args.cloud_product)
        if cloud_data.data_parity==False: continue
        trop_data.calc_geo_column()
        trop_data.apply_bias_correction()
        trop_data.cloud_filter_and_preprocess(cloud_data, cloud_threshold, args.pmax, args.pmin)
        grid_aggregator.initalise_grid()
        grid_aggregator.grid_trop_data(trop_data)
        grid_aggregator.apply_cloud_slice()
    grid_aggregator.calc_seasonal_means()

    out_file = path.join(args.out_dir, 'tropomi-ut-no2-'+args.cloud_product
                         + '-' + args.cloud_threshold
                         + '-' + args.grid_res
                         + '-' + args.season
                         + '-' + yrrange+'-v2.nc')
    grid_aggregator.print_report()
    grid_aggregator.save_to_netcdf(out_file)
    grid_aggregator.plot_data()




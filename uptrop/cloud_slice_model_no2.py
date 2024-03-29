#!/usr/bin/python

'''
Use synthetic partial columns from GEOS-Chem to obtain cloud-sliced
NO2 in the upper troposphere and compare this to UT NO2 obtained if
simply average the NO2 mixing ratios from the model over the same
pressure range (the "truth"). Both are obtained by Gaussian
weighting toward the pressure center.

GEOS-Chem partial columns are obtained over Europe, North America,
and China at the GEOS-FP meteorology native resolution (0.25x0.3125)
(latxlon) for June-August 2016-2017.

Input options to process the data include the region, the horizontal
resolution, and model simulation years.

.. code-block:: bash

    usage: ut_no2_gc_test.py [-h] [--gc_dir GC_DIR] [--out_dir OUT_DIR]
                         [--resolution RESOLUTION] [--region REGION]
                         [--strat_filter_threshold STRAT_FILTER_THRESHOLD]
                         [--start_date START_DATE]
                         [--end_date END_DATE] [-p PLOT]
                         [--do_temp_correct DO_TEMP_CORRECT]
                         [--apply_cld_frac_filter APPLY_CLD_FRAC_FILTER]
                         [--do_cld_hght_test DO_CLD_HGHT_TEST]

    optional arguments:
      -h, --help            show this help message and exit
      --gc_dir GC_DIR
      --out_dir OUT_DIR
      --resolution RESOLUTION
                            Can be 8x10, 4x5, 2x25 or 1x1
      --region REGION       Can be EU, NA, or CH
      --strat_filter_threshold STRAT_FILTER_THRESHOLD
      --start_date START_DATE
      --end_date END_DATE
      -p PLOT, --plot PLOT
      --do_temp_correct DO_TEMP_CORRECT
      --apply_cld_frac_filter APPLY_CLD_FRAC_FILTER
      --do_cld_hght_test DO_CLD_HGHT_TEST

'''

# Note for the future; most of this can probably be replaced with a modification of the GridAggregator class from
# tropomi_ut_no2

# Import relevant packages:
import numpy as np
from netCDF4 import Dataset
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import argparse
from sklearn.linear_model import LinearRegression
import sys
import os
import datetime as dt
from dateutil import rrule as rr

# Import hack
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..'))

from uptrop.date_file_utils import get_gc_file_list
from uptrop.constants import AVOGADRO
from uptrop.constants import G
from uptrop.constants import MW_AIR
from uptrop.bootstrap import rma
from uptrop.cloud_slice_no2 import cldslice, CLOUD_SLICE_ERROR_ENUM
from uptrop.height_pressure_converter import alt2pres, pres2alt


# Turn off warnings:
np.warnings.filterwarnings('ignore')


# Define pressure range:
# These are fixed, so can be constants, as opposed to inputs.
P_MIN=180     
P_MAX=450     


class ProcessingException(Exception):
    pass


class CloudSliceException(Exception):
    pass


class InvalidRegionException(Exception):
    pass


class InvalidResolutionException(Exception):
    pass


class DomainIssueException(Exception):
    pass


class ProcessedData:
    """A class for comparing geoschem and tropomi data"""
    # Note for anyone reading this; this performs very similar functions to the GridAggregator class, but in a different
    # way. I prefer the method in GridAggregator, but this is fine for now. -John Roberts
    def __init__(self, region, str_res, strat_thld, do_temperature_correction=False, do_cld_frac_filter=False, do_cld_hght_test=False):
        """Creates and returns an instance of ProcessedData for a given region

        This class aggregates geoschem data over a grid defined by the :ref:`define_grid method`.

        :param region: The region to aggregate geoschem over; can be NA, EU or CH.
        :type region: str
        :param str_res: The resolution of the grid around teh region to analyse; can be 8x10, 4x5, 2x25 or 1x1
        :type str_res: str
        :param do_temperature_correction: Whether to perform the temperature correction step
        :type do_temperature_correction: bool
        :param do_cld_frac_filter: Whether to perform fractional filtering of clouds
        :type do_cld_frac_filter: bool

        :returns: A ProcessedData isntance ready to recieve data
        :rtype: ProcessedData
        """

        self.temperature_correction = do_temperature_correction
        self.cld_frac_filter = do_cld_frac_filter
        self.cloud_height_test = do_cld_hght_test
        
        # Filtering threshold for stratosphere:
        self.strat_filt = strat_thld

        self.define_grid(region, str_res)
        grid_shape = self.X.shape

        # Define output arrays:
        self.g_no2_vcd = np.zeros(grid_shape)
        self.g_no2_vmr = np.zeros(grid_shape)
        self.g_cld_fr = np.zeros(grid_shape)
        self.g_cld_p = np.zeros(grid_shape)
        self.g_slope_err = np.zeros(grid_shape)
        self.g_gaus_wgt = np.zeros(grid_shape)
        self.true_o3 = np.zeros(grid_shape)  # ozone mixing ratio coincident with cloud-sliced ut self
        self.true_no2 = np.zeros(grid_shape)  # "true" cloudy UT NO2
        self.g_askut_no2 = np.zeros(grid_shape)  # "true" all-sky UT NO2
        self.g_ask_gaus_wgt = np.zeros(grid_shape)  # "true" all-sky UT NO2 weights
        self.g_as_cnt = np.zeros(grid_shape)  # Count all-sky
        self.g_cnt = np.zeros(grid_shape)

        # Type of data loss in the cloud-slicing retrieval:
        self.loss_count = {
            "too_few_points": 0,
            "low_cloud_height_range": 0,
            "low_cloud_height_std": 0,
            "large_error": 0,
            "sig_diff_from_zero": 0,
            "no2_outlier": 0,
            "non_uni_strat": 0,
        }

        # Initialize:
        self.cloud_slice_count = 0
        self.maxcnt = 0
        self.grad_retain = 0
        self.grad_remove = 0
        # Print out min and max cloud-sliced UT NO2 error size:
        self.minerr = 1   # Set to maximum possible error (100%)
        self.maxerr = 0

        # Define string to represent the layer range:
        self.prange = str(P_MIN) + '-' + str(P_MAX)

        # Define factor to convert slope of NO2 mixing ratio versus pressure
        # to VMR:
        self.den2mr = np.divide((np.multiply(G, MW_AIR)), AVOGADRO)

    def define_grid(self, region, str_res):
        """Defines a grid based on a set of regional bounding boxes and a resolution

        :param region: The region to aggregate geoschem over; can be NA, EU or CH.
        :type region: str
        :param str_res: The resolution of the grid around the region to analyse; can be 8x10, 4x5, 2x25 or 1x1
        :type str_res: str
        """
        # Define target grid:
        if region == 'NA':
            self.minlat = 2.
            self.maxlat = 62.
            self.minlon = -140.
            self.maxlon = -55.
        elif region == 'EU':
            self.minlat = 26.
            self.maxlat = 66.
            self.minlon = -25.
            self.maxlon = 45.
        elif region == 'CH':
            self.minlat = 6.
            self.maxlat = 62.
            self.minlon = 60.
            self.maxlon = 140.
        else:
            print("Invalid region; valid regions are 'NA','EU','CH'.")
            raise InvalidRegionException
        # NOTE FOR LATER: This snippet turns up in fresco_cld_err; a candidate for the library.
        # Define grid information:
        if str_res == '8x10':
            self.dellat, self.dellon = 8, 10
        elif str_res == '4x5':
            self.dellat, self.dellon = 4, 5
        elif str_res == '2x25':
            self.dellat, self.dellon = 2, 2.5
        elif str_res == '1x1':
            self.dellat, self.dellon = 1, 1
        else:
            print("Invalid resolution: valid resolutions are 8x10, 4x5, 2x25 (two pt 5) and 1x1")
            raise InvalidResolutionException
        self.out_lon = np.arange(self.minlon, self.maxlon + self.dellon, self.dellon)
        self.out_lat = np.arange(self.minlat, self.maxlat + self.dellat, self.dellat)
        # Convert output lats and long to 2D:
        self.X, self.Y = np.meshgrid(self.out_lon, self.out_lat, indexing='ij')
        # Dimensions of output data:
        self.xdim = len(self.out_lon)
        self.ydim = len(self.out_lat)
        # Get maximum possible number of native resolution model gridsquares
        # in the coarser grid (to use in process_grid_square to check that 
        # don't exceed this):
        self.max_limit=(self.dellon/0.3125)*(self.dellat/0.25)


    # ----Processing methods----
    def process_geoschem_day(self, file_path):
        """Aggregates geoschem data into the processing grid for a file on a given day.

        Calls :ref:`regrid_and_process` for each pixel in the geoschem file, then calls
        :ref:`process_grid_square` for each square of the grid

        :param file_path: Path to the geoschem file to load into the processing grid
        :type file_path: str"""
        # Define output data for this day:
        out_shape = (self.xdim, self.ydim)  # This feel gross. A 3-d list of appendable lists.
        self.g_no2 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.g_grad = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.g_o3 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.all_cld_fr = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.strat_no2 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.all_cld_p = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.g_true_no2 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]

        this_geoschem_day = GeosChemDay(file_path,
                                        temperature_correction=self.temperature_correction,
                                        cloud_height_test=self.cloud_height_test)
        # Get column values:
        for y in range(len(this_geoschem_day.t_lat)):
            for x in range(len(this_geoschem_day.t_lon)):
                this_geoschem_day.prepare_no2_pixel(x, y)
                self.regrid_and_process(x, y, this_geoschem_day)
        for i in range(self.xdim):
            for j in range(self.ydim):
                if any(value != 0 for value in self.g_no2[i][j]):
                    self.process_grid_square(i, j)

    def regrid_and_process(self, x, y, gc_data):
        """Inserts data from pixel x,y from geoschem data into the appropriate grid cell

        :param x: X index of geoschem grid cell
        :type x: int
        :param y: Y index of geoschem grid cell
        :type y: int
        :param gc_data: An instance of GeosChemDay data
        :type gc_data: GeosChemDay
        """
        # If no valid data, skip:
        if ( len(gc_data.askind) == 0):
            return

        # Find nearest gridsquare in output grid:
        lon = int(np.argmin(abs(self.out_lon - gc_data.t_lon[x])))
        lat = int(np.argmin(abs(self.out_lat - gc_data.t_lat[y])))

        self.g_ask_gaus_wgt[lon, lat] += np.sum(gc_data.twgt)

        self.g_askut_no2[lon, lat] += np.sum(gc_data.t_gc_no2[gc_data.askind, y, x] * gc_data.twgt * 1e3)
        self.g_as_cnt[lon, lat] += 1.0

        # Add relevant data to array if there are clouds at 450-180 hPa:
        if (gc_data.lcld < gc_data.level_min) or (gc_data.lcld > gc_data.level_max):
            # Check for clouds between min and max pressure. If none, move to next pixel.
            #print("Cloud top outside pressure range in pixel {},{}".format(x,y))
            return

        # Find where cloud fraction in UT exceeds 0.7 after calculating
        # true all-sky NO2:
        # (Keep for testing effect of thick clouds on cloud-sliced UT NO2):
        if ( self.cld_frac_filter ):
            if np.sum(gc_data.t_cld_fr[gc_data.level_min:gc_data.level_max + 1, y, x]) < 0.7:
                return

        self.g_no2[lon][lat].append(gc_data.no2_2d)
        self.g_grad[lon][lat].append(gc_data.grad)
        self.strat_no2[lon][lat].append(gc_data.strat_col)
        self.all_cld_p[lon][lat].append(gc_data.t_cld_hgt[y, x])
        self.g_o3[lon][lat].append(
            np.mean(gc_data.t_gc_o3[gc_data.level_min:gc_data.level_max + 1, y, x]))
        self.g_true_no2[lon][lat].append(
            np.mean(gc_data.t_gc_no2[gc_data.level_min:gc_data.level_max + 1, y, x]))
        self.all_cld_fr[lon][lat].append(
            np.sum(gc_data.t_cld_fr[gc_data.level_min:gc_data.level_max + 1, y, x]))
        pass

    def process_grid_square(self, i, j):

        # Define vectors of relevant data:
        # These should all have the same length
        t_col_no2 = np.asarray(self.g_no2[i][j],dtype=np.float)
        t_strat_no2 = np.asarray(self.strat_no2[i][j],dtype=np.float)
        t_fr_c = np.asarray(self.all_cld_fr[i][j],dtype=np.float)
        t_cld = np.asarray(self.all_cld_p[i][j],dtype=np.float)
        t_mr_no2 = np.asarray(self.g_true_no2[i][j],dtype=np.float)
        t_o3 = np.asarray(self.g_o3[i][j],dtype=np.float)
        t_grad_no2 = np.asarray(self.g_grad[i][j],dtype=np.float)

        # Skip if fewer than 10 points:
        if len(t_col_no2) < 10:
            #print("Grid square {}, {} failed with condition too_few_points".format(i, j))
            self.loss_count["too_few_points"] += 1
            return

        # Remove non-uniform stratosphere:
        if (np.std(t_strat_no2) / np.mean(t_strat_no2)) > self.strat_filt:
            self.loss_count["non_uni_strat"] += 1
            return  # continue

        # Get number of points:
        n_pnts = len(t_cld)
        if n_pnts > self.maxcnt:
            self.maxcnt = n_pnts
            print('Max no. of points in a cluster: ', self.maxcnt, flush=True)

        # Check if number of gridsquares exceeds max possible
        # (indicates error in coarser grid domain):
        if ( self.maxcnt > self.max_limit ):
            print('No. of model grids exceeds max possible')
            raise DomainIssueException

        # Use cloud_slice_ut_no2 function to get cloud-sliced UT NO2.
        # Treat data differently depending on whether there are 10-99 points:
        if n_pnts >= 10 and n_pnts <100:
            self.add_slice(i,j,t_cld,t_col_no2, t_mr_no2, t_fr_c, t_grad_no2, t_o3)
        # Or at least 100 points:
        elif n_pnts >= 100:
            num_slices = 40
            stride = round(n_pnts / num_slices)
            nloop = list(range(stride))
            for w in nloop:
                subset_t_col_no2 = t_col_no2[w::stride]
                subset_t_cld = t_cld[w::stride]
                subset_t_mr_no2 = t_mr_no2[w::stride]
                subset_t_fr_c = t_fr_c[w::stride]
                subset_t_grad_no2 = t_grad_no2[w::stride]
                subset_t_o3 = t_o3[w::stride]
                self.add_slice(i, j, subset_t_cld, subset_t_col_no2, subset_t_mr_no2, subset_t_fr_c, subset_t_grad_no2, subset_t_o3)

    def add_slice(self, i, j, t_cld, t_col_no2, t_mr_no2, t_fr_c, t_grad_no2, t_o3):
        """Extracts the upper troposphere gc_data, gc_data error, ozone data, ozone error,
         and mean cloud pressure for grid square [i,j]

        This method uses the cloud-slicing function :ref:`uptrop.cloud_slice_ut_no2.cldslice`
        Once calculated, the a weighting is derived from cloud pressure.
        The weighted upper tropospheric gc_data and error is added to the rolling total for this season.
        If the cloud slicing fails, then the reason is added to loss_count for the end report.

        :param i: X-index of grid square
        :type i: int
        :param j: Y-index of grid square
        :type j: int
        :param t_cld: Tropospheric cloud
        :type t_cld: float
        :param t_col_no2: Total (stratosphere+troposphere) column no2
        :type t_col_no2: float
        :param t_mr_no2: NO2 mixing ratio between pmax and pmin
        :type t_mr_no2: float
        :param t_fr_c: Cloud fraction between pmax and pmin
        :type t_fr_c: float
        :param t_grad_no2: NO2 gradient between pmax and pmin in pptv/hPa
        :type t_grad_no2: float
        :param t_o3: O3 mixing ratio to test influence from stratosphere
        :type t_o3: float
        """
        utmrno2, utmrno2err, stage_reached, mean_cld_pres = cldslice(t_col_no2, t_cld,140)
        
        # Skip if approach didn't work (i.e. cloud-sliced UT NO2 is NaN):
        # Drop out after the reason for data loss is added to loss_count.
        if np.isnan(utmrno2) or np.isnan(utmrno2err):
            # Cloud-slicing led to nan, but due to rma regression rather
            # than data filtering (these are rare):
            if ( stage_reached==0 ):
                print("Cloud-sliced NO2 NAN for pixel i:{} j:{}".format(i, j))
                return
            self.loss_count[CLOUD_SLICE_ERROR_ENUM[stage_reached]] += 1
            # Track non-uniform NO2 removed:
            grad_ind=np.where(np.abs(t_grad_no2) >= 0.33)[0]
            self.grad_remove += len(grad_ind)
            #print("Cloud-slice exception {} in pixel i:{} j:{}".format(
            #    CLOUD_SLICE_ERROR_ENUM[stage_reached], i, j))
        else:
            # Calculate weights:
            #err_wgt = 1.0 / (utmrno2err ** 2)  # not use, but preserve anyway
            gaus_wgt = np.exp((-(mean_cld_pres - 315) ** 2) / (2 * 135 ** 2))
            
            # Get and track error range:
            nerr = utmrno2err / utmrno2
            if nerr > self.maxerr:
                self.maxerr = nerr
                print('Max rel. cloud-sliced NO2 error: ', self.maxerr, flush=True)
            if nerr < self.minerr:
                self.minerr = nerr
                print('Min rel. cloud-sliced NO2 error: ', self.minerr, flush=True)
                
            # Track non-uniform NO2 retained:
            grad_ind=np.where(np.abs(t_grad_no2) >= 0.33)[0]
            self.grad_retain += len(grad_ind)
            
            # Weighted mean for each pass of cldslice:
            weight_mean = np.mean(t_mr_no2 * 1e3) * gaus_wgt 
            self.true_no2[i, j] += weight_mean
            self.true_o3[i, j] += (np.mean(t_o3) * gaus_wgt)
            self.g_no2_vmr[i, j] += utmrno2 * gaus_wgt
            self.g_slope_err[i, j] += utmrno2err * gaus_wgt
            self.g_gaus_wgt[i, j] += gaus_wgt
            self.g_cnt[i, j] += 1
            self.g_cld_fr[i, j] += np.mean(t_fr_c)
            self.cloud_slice_count += 1

    def get_weighted_mean(self):
        """Applies weighting to the aggregated means.
        """
        # Mean physical parameters:
        self.g_no2_vmr = np.divide(self.g_no2_vmr, self.g_gaus_wgt, where=self.g_cnt != 0)
        self.g_cld_fr = np.divide(self.g_cld_fr, self.g_cnt, where=self.g_cnt != 0)
        self.true_no2 = np.divide(self.true_no2, self.g_gaus_wgt, where=self.g_cnt != 0)
        self.g_cld_p = np.divide(self.g_cld_p, self.g_cnt, where=self.g_cnt != 0)
        self.true_o3 = np.divide(self.true_o3, self.g_gaus_wgt, where=self.g_cnt != 0)
        self.g_askut_no2 = np.divide(self.g_askut_no2, self.g_ask_gaus_wgt, where=self.g_as_cnt != 0)
        # Mean cloud-sliced UT NO2 error:
        self.g_slope_err = np.divide(self.g_slope_err, self.g_gaus_wgt, where=self.g_cnt != 0)
        # Mean weights:
        self.g_gaus_wgt = np.divide(self.g_gaus_wgt, self.g_cnt, where=self.g_cnt != 0)
        # No data (nan):
        self.true_no2[self.g_cnt == 0] = np.nan
        self.true_o3[self.g_cnt == 0] = np.nan
        self.g_slope_err[self.g_cnt == 0] = np.nan
        self.g_gaus_wgt[self.g_cnt == 0] = np.nan
        self.g_askut_no2[self.g_as_cnt == 0] = np.nan
        self.g_no2_vmr[self.g_cnt == 0] = np.nan
        self.g_cld_fr[self.g_cnt == 0] = np.nan
        self.g_cld_p[self.g_cnt == 0] = np.nan

    # ----Reporting and saving methods----
    def print_data_report(self):
        """Prints a set of reasons for missing data
        """
        # Output code diagnostics:
        # No. of data points:
        # print('No. of valid data points: ',cloud_slice_count,flush=True)
        print('Max no. of data points in a gridsquare: ', np.amax(self.g_cnt), flush=True)
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
        print('Non-uniform NO2 retained = ', self.grad_retain)
        print('Non-uniform NO2 removed = ', self.grad_remove)
        print('Error range of cloud-sliced UT NO2: ', np.nanmin(self.g_slope_err), np.nanmax(self.g_slope_err))

    def plot_data(self):
        """Plots the present state of the gridded data
        """
        # ===> FUTURE UPDATE : change from basemap to cartopy <===
        # Plot the data:
        m = Basemap(resolution='l', projection='merc',
                    lat_0=0, lon_0=0, llcrnrlon=self.minlon,
                    llcrnrlat=self.minlat, urcrnrlon=self.maxlon, urcrnrlat=self.maxlat)
        xi, yi = m(self.X, self.Y)
        plt.subplot(2, 3, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.g_no2_vmr), vmin=0, vmax=80, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('Cloud-sliced NO2 VMRs')
        plt.subplot(2, 3, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.true_no2), vmin=0.0, vmax=80, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('True cloudy NO2')
        plt.subplot(2, 3, 3)
        cs = m.pcolor(xi, yi, np.squeeze(self.g_askut_no2), vmin=0, vmax=80, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('True NO2 VMRs under all-sky conditions')
        plt.subplot(2, 3, 4)
        plt.plot(self.true_no2, self.g_no2_vmr, 'o', color='black', markersize=6)
        r = stats.pearsonr(self.true_no2[~np.isnan(self.g_no2_vmr)],
                           self.g_no2_vmr[~np.isnan(self.g_no2_vmr)])
        # print('Correlation = ', r[0])
        result = rma(self.true_no2[~np.isnan(self.g_no2_vmr)], self.g_no2_vmr[~np.isnan(self.g_no2_vmr)],
                     len(self.true_no2[~np.isnan(self.g_no2_vmr)]), 1000)
        print(result, flush=True)
        xvals = np.arange(0, 100, 5)
        yvals = result[1] + xvals * result[0]
        plt.plot(xvals, yvals, '-')
        plt.xlim(-4, 80)
        plt.ylim(-4, 80)
        plt.xlabel('True NO2 (cloudy)')
        plt.ylabel('Cloud-sliced NO2')
        print('===== True (cloudy) vs cloud-sliced UT NO2 ====')
        print('R = ', r[0], flush=True)
        print('Slope = ', result[0])
        print('Slope Err = ', result[2], flush=True)
        print('Intercept = ', result[1], flush=True)
        print('Intercept Err = ', result[3], flush=True)
        add2plt = ("y = {a:6.2f}x + {b:6.3f}".format(a=result[0], b=result[1]))
        plt.text(2, 75, add2plt, fontsize=8,
                 ha='left', va='center')  # , transform=ax.transAxes)
        add2plt = ("R = {a:6.2f}".format(a=r[0]))
        plt.text(2, 65, add2plt, fontsize=8,
                 ha='left', va='center')  # , transform=ax.transAxes)
        plt.subplot(2, 3, 5)
        plt.plot(self.g_askut_no2, self.g_no2_vmr, 'o', color='black', markersize=6)
        r = stats.pearsonr(self.g_askut_no2[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_askut_no2))],
                           self.g_no2_vmr[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_askut_no2))])
        result = rma(self.g_askut_no2[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_askut_no2))],
                     self.g_no2_vmr[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_askut_no2))],
                     len(self.g_askut_no2[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_askut_no2))]),
                     1000)
        xvals = np.arange(0, 100, 5)
        yvals = result[1] + xvals * result[0]
        plt.plot(xvals, yvals, '-')
        plt.xlim(-4, 80)
        plt.ylim(-4, 80)
        plt.xlabel('True NO2 (all-sky)')
        plt.ylabel('Cloud-sliced NO2')
        add2plt = ("y = {a:6.2f}x + {b:6.3f}". \
                   format(a=result[0], b=result[1]))
        plt.text(2, 75, add2plt, fontsize=8, \
                 ha='left', va='center')  # , transform=ax.transAxes)
        add2plt = ("R = {a:6.2f}".format(a=r[0]))
        plt.text(2, 65, add2plt, fontsize=8, \
                 ha='left', va='center')  # , transform=ax.transAxes)
        plt.subplot(2, 3, 6)
        plt.plot(self.g_askut_no2, self.true_no2, 'o', color='black', markersize=6)
        r = stats.pearsonr(self.g_askut_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_askut_no2))],
                           self.true_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_askut_no2))])
        result = rma(self.g_askut_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_askut_no2))],
                     self.true_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_askut_no2))],
                     len(self.g_askut_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_askut_no2))]),
                     1000)
        xvals = np.arange(0, 100, 5)
        yvals = result[1] + xvals * result[0]
        plt.plot(xvals, yvals, '-')
        plt.xlim(-4, 80)
        plt.ylim(-4, 80)
        plt.xlabel('True NO2 (all-sky)')
        plt.ylabel('True NO2 (cloudy)')
        add2plt = ("y = {a:6.2f}x + {b:6.3f}".format(a=result[0], b=result[1]))
        plt.text(2, 75, add2plt, fontsize=8, ha='left', va='center')
        add2plt = ("R = {a:6.2f}".format(a=r[0]))
        plt.text(2, 65, add2plt, fontsize=8, ha='left', va='center')
        plt.show()

    def save_to_netcdf(self, out_path):
        """Saves the present state of the grid to a netCDF4 file

        :param out_path: Path to the output file
        :type out_path: str
        """
        # Save the data to NetCDF:
        ncout = Dataset(out_path, mode='w', format='NETCDF4')
        # Create data dimensions:
        ncout.createDimension('lat', self.ydim)
        ncout.createDimension('lon', self.xdim)
        # create lon axis:
        lon = ncout.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        lon[:] = self.out_lon
        # create lat axis:
        lat = ncout.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latgitude'
        lat[:] = self.out_lat
        # create data axes:
        # (1) Cloud-sliced NO2:
        csutno2 = ncout.createVariable('csutno2', np.float32, ('lon', 'lat'))
        csutno2.units = 'pptv'
        csutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained using cloud-slicing'
        csutno2[:] = self.g_no2_vmr
        # (2a) Double-weighting error:
        utdblerr = ncout.createVariable('utdblerr', np.float32, ('lon', 'lat'))
        utdblerr.units = 'pptv'
        utdblerr.long_name = 'Standard error of the NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
        utdblerr[:] = self.g_slope_err
        # (2b) Gaussian-weighting error:
        utgauserr = ncout.createVariable('utgauserr', np.float32, ('lon', 'lat'))
        utgauserr.units = 'pptv'
        utgauserr.long_name = 'Standard error of the NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
        utgauserr[:] = self.g_gaus_wgt
        # (3) Number of observations in each gridsquare:
        nobs = ncout.createVariable('nobs', np.float32, ('lon', 'lat'))
        nobs.units = 'unitless'
        nobs.long_name = 'Number of observations in each gridsquare used to obtain cloud-sliced UT NO2 mixing ratios'
        nobs[:] = self.g_cnt
        # (4) Mean cloud pressure for season between 450-180 hPa:
        utcld = ncout.createVariable('utcld', np.float32, ('lon', 'lat'))
        utcld.units = 'hPa'
        utcld.long_name = 'Mean cloud pressure between 450 and 180 hPa'
        utcld[:] = self.g_cld_p
        # (5) Mean NO2 mixing ratio at 450-180 hPa for scenes with clouds:
        cldutno2 = ncout.createVariable('cldutno2', np.float32, ('lon', 'lat'))
        cldutno2.units = 'pptv'
        cldutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained if clouds are present'
        cldutno2[:] = self.true_no2
        # (6) Mean NO2 mixing ratio at 450-180 hPa under all conditions (all-sky):
        askutno2 = ncout.createVariable('askutno2', np.float32, ('lon', 'lat'))
        askutno2.units = 'pptv'
        askutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained under all conditions (all-sky)'
        askutno2[:] = self.g_askut_no2
        # (7) Cloud fraction:
        utcldfrc = ncout.createVariable('utcldfrc', np.float32, ('lon', 'lat'))
        utcldfrc.units = 'unitless'
        utcldfrc.long_name = 'GEOS-FP cloud fraction obtained as sum of 3D cloud fractions across range of interest (180-450 hPa)'
        utcldfrc[:] = self.g_cld_fr
        # (8) O3 sampled coincident with cloud-slicing retrieval:
        uto3 = ncout.createVariable('uto3', np.float32, ('lon', 'lat'))
        uto3.units = 'ppbv'
        uto3.long_name = 'GEOS-Chem ozone obtained coincident with cloud-sliced NO2'
        uto3[:] = self.true_o3
        # Close the file:
        ncout.close()


class GeosChemDay:
    """A class for reading, preprocessing and accessing Geoschem data on a given day
    """
    def __init__(self, file_path, temperature_correction=False, cloud_height_test=False):
        """Reads the data at file_path and returns a GeosChemDay object containing that data

        :param file_path: Path to the netcdf4 file containing the GeosChem data
        :type file_path: str
        :param temperature_correction: Whether to apply temperature correction
        :type temperature_correction: bool
        :param cloud_height_test: Whether to test effect of systematic underestimate in cloud height
        :type cloud_height_test: bool

        :returns: A GeosChemDay class
        :rtype: GeosChemDay
        """
        print('File path: ',file_path, flush=True)

        self.temperature_correction = temperature_correction
        self.cloud_height_test = cloud_height_test

        # Read dataset:
        fh = Dataset(file_path, mode='r')
        # Extract data of interest:
        tlon, tlat, tgcno2, tcldfr, tcldhgt, tadn, tbxhgt, tpedge, tpause, tgco3, tdegk = \
            fh.variables['LON'], fh.variables['LAT'], \
            fh.variables['IJ-AVG-S__NO2'], fh.variables['TIME-SER__CF'], \
            fh.variables['TIME-SER__CThgt'], fh.variables['TIME-SER__AIRDEN'], \
            fh.variables['BXHGHT-S__BXHEIGHT'], fh.variables['PEDGE-S__PSURF'], \
            fh.variables['TR-PAUSE__TP-PRESS'], fh.variables['IJ-AVG-S__O3'], \
            fh.variables['DAO-3D-S__TMPU']
        self.t_lon = tlon[:]
        self.t_lat = tlat[:]
        self.t_gc_no2 = tgcno2[:]
        self.t_cld_fr = tcldfr[:]
        self.t_cld_hgt = tcldhgt[0, :, :]
        self.t_adn = tadn[:]  # in molec/cm3
        self.t_bx_hgt = tbxhgt[:]
        self.t_p_edge = tpedge[:]
        self.t_pause = tpause[0, :, :]
        self.t_gc_o3 = tgco3[:]
        self.t_deg_k = tdegk[:]
        # Convert box height from m to cm:
        self.t_bx_hgt = self.t_bx_hgt * 1e2

        if self.cloud_height_test:
            # Lower cloud heights by 1 km to roughly mimic lower altitude
            # clouds retrieved for TROPOMI assuming clouds are reflective
            # boundaries with uniform reflectivity:
            # Calculate cloud top height in m:
            t_cld_hgt = pres2alt(self.t_cld_hgt*1e2)
            # Lower the clouds by 1 km (1000 m) (this won't work for low-altitude
            # clouds):
            t_cld_hgt = t_cld_hgt - 1e3
            # Convert back to Pa and convert that to hPa:
            self.t_cld_hgt = alt2pres(t_cld_hgt)*1e-2

        #Get outputs ready here for tidyness:
        self.no2_2d = None
        self.trop_col = None
        self.strat_col = None
        self.gcutno2 = None
        self.gascnt = None
        self.grad = None

        self.level_min = None
        self.level_max = None
        self.askind = None

    def prepare_no2_pixel(self, x, y):
        """Extracts preprocesed no2 from the geoschem pixel at x,y

        :param x: The x index of the pixel
        :type x: int
        :param y: The y index of the pixel
        :type y: int
        """
        # Calculate corresponding mid-pressure values:
        tp_mid = np.zeros(len(self.t_p_edge[:, y, x]))
        # Get mid-pressure values, except for highest layer:
        for k in range(len(self.t_p_edge[:, y, x]) - 1):
            tp_mid[k] = np.multiply(0.5, (self.t_p_edge[k, y, x] + self.t_p_edge[k + 1, y, x]))
        # Estimate mid-pressure for highest value (doesn't need to
        # be accurate, as surpassing the range of interset):
        # Data output from the model includes 47 vertical layers. This means that only 46 pressure centres can be calculated as the calculation requires pressure edges.
        tp_mid[46] = np.multiply(0.5, (self.t_p_edge[46, y, x] + (self.t_p_edge[46, y, x] - 0.1)))
        # Get model layer of tropopause:
        tppind = np.argmin(abs(tp_mid - self.t_pause[y, x]))
        # Get indices that fall between 450 and 180 hPa for estimating
        # "true' all-sky UT NO2 and partial columns:
        lind = np.where((tp_mid >= P_MIN) & (tp_mid <= P_MAX))[0]
        # Get UT NO2 under "true" all-sky conditions:
        # Make sure this is below the tropopause:
        # If below tropopause, use full extent (180-450 hPa):
        if lind[len(lind) - 1] <= tppind:
            self.askind = lind
        # If above tropopause, trim to tropopause-450 hPa:
        if lind[len(lind) - 1] > tppind:
            self.askind = lind[np.where(lind <= tppind)[0]]
        # If tropopause below 450 hPa, skip entirely:
        if self.t_pause[y, x] > P_MAX:
            #print("Tropopause less than P_MAX in geoschem pixel x:{}, y:{}".format(x,y))
            return  # continue

        # Get Guassian weights that allocate higher weights to points
        # closest to the pressure centre (315 hPa):
        # Equation is:
        #   w = exp(-(p-315)^2/2*135^2 ) where 315 hPa is the centre and
        #         135 hPa is the standard deviation.
        self.twgt = np.exp((-(tp_mid[self.askind] - 315) ** 2) / (2 * 135 ** 2))

        # Get model level of cloud top height closest to lowest
        # pressure-altitude of interest (P_MIN):
        self.lcld = np.argmin(abs(self.t_cld_hgt[y, x] - tp_mid))
        # Skip if cloud top height ouside pressure range of interest:
        self.level_min, self.level_max = np.amin(lind), np.amax(lind)

        if (self.temperature_correction):

            # Equation is from the TROPOMI product ATBD (p. 32, Eqn 18)
            # (product document abbrevation: S5P-KNMI-L2-0005-RP)
            self.temp_corr = 1 - (3.16e-3 * (self.t_deg_k[self.level_min:, y, x] - 220.)) + \
                        (3.39e-6 * ((self.t_deg_k[self.level_min:, y, x] - 220) ** 2))
        else:
            # Set to 1 so that no scaling is applied:
            # (might be a more eloquent way to do this)
            self.temp_corr = np.ones(len(self.t_gc_no2[self.level_min:, y, x]))

        # Calculate NO2 gradient:
        regr = LinearRegression()  
        # Pressure (hPa)
        x_vals=tp_mid[self.level_min:self.level_max].reshape(-1,1)
        # NO2 (ppbv)
        y_vals=self.t_gc_no2[self.level_min:self.level_max, y, x].reshape(-1,1)
        # NO2 (pptv)
        y_vals=y_vals*1e3
        # Perform regression:
        regr.fit(x_vals,y_vals)
        # Define gradient from regression slope (pptv/hPa):
        self.grad = regr.coef_

        # Get partial NO2 column in molec/m2 from cloud top height
        # to highest model level (output up to level 47):
        # print(t_gc_no2[self.level_min:tppind,y,x])
        # print(t_gc_no2[self.level_min:tppind,y,x]*1.5)
        self.no2_2d = np.sum(self.t_gc_no2[self.level_min:, y, x]
                             * 1e-5
                             * self.temp_corr
                             * self.t_adn[self.level_min:, y, x]
                             * self.t_bx_hgt[self.level_min:, y, x])
        # Get stratospheric column from 180 hPa aloft:
        # Previous approach (remove when model simulations done):
        # tppind=np.where(tpmid<180.)[0]
        self.strat_col = np.sum(self.t_gc_no2[tppind:, y, x]
                                * 1e-5 * self.t_adn[tppind:, y, x]
                                * self.t_bx_hgt[tppind:, y, x])


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    # Shorten directory name to up to "GC/", then define the subdirectory
    # as 'geosfp' + dirreg + 'iccw/' in get_file_list.
    # This is now done in get_gc_file_list
    parser.add_argument("--gc_dir")
    parser.add_argument("--out_dir")
    parser.add_argument("--resolution", default="4x5", help="Can be 8x10, 4x5, 2x25 or 1x1")
    parser.add_argument("--region", default="EU", help="Can be EU, NA, or CH")
    parser.add_argument("--strat_filter_threshold", default="002", help="")
    #parser.add_argument("--start_date", default="2016-06-01")
    #parser.add_argument("--end_date", default="2017-08-31")
    parser.add_argument("-p", "--plot", type=bool, default=False)
    parser.add_argument("--do_temp_correct", type=bool)
    parser.add_argument("--apply_cld_frac_filter", type=bool)
    parser.add_argument("--do_cld_hght_test", type=bool)
    args = parser.parse_args()

    # Get files:
    gc_dir = args.gc_dir
    STR_RES = args.resolution
    REGION = args.region
    
    # Hard code for GEOS-Chem synthetic experiment:
    yrrange = '2016-2017'

    # Define output file depending on input arguments:
    out_file = os.path.join(args.out_dir, 'gc-v12-1-0-ut-no2'
                            + '-' + args.region
                            + '-' + args.resolution
                            + '-jja-' + yrrange
                            + '-' + "gaus-wgt" )

    if args.do_cld_hght_test:
        out_file += "-cldtop"
    if args.strat_filter_threshold != "002":
        out_file += '-' + args.strat_filter_threshold +'strat'
    if args.do_temp_correct:
        out_file += "-temp-corr"
    if args.apply_cld_frac_filter:
        out_file += "-cld-filt"
        
    out_file += "-v1.nc"

    strat_filter_threshold = float(args.strat_filter_threshold)/100
    files = get_gc_file_list(gc_dir, args.region)
    print('Number of files:', len(files), flush=True)

    rolling_total = ProcessedData(REGION, STR_RES, strat_filter_threshold,
                                  do_temperature_correction=args.do_temp_correct,
                                  do_cld_frac_filter=args.apply_cld_frac_filter,
                                  do_cld_hght_test=args.do_cld_hght_test)

    for file_path in files:
        rolling_total.process_geoschem_day(file_path)

    rolling_total.get_weighted_mean()
    rolling_total.print_data_report()
    if args.plot==True: rolling_total.plot_data()
    rolling_total.save_to_netcdf(out_file)



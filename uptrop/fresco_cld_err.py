#!/usr/bin/python

''' Read and regrid cloud information from the official TROPOMI CLOUD product 
    and from the FRESCO cloud product in the TROPOMI NO2 data product file.

    Output from this code includes a NetCDF file with:
    (1) Data regridded to a 1 deg x 1 deg global grid and saved as monthly 
        means for scenes with the reference cloud product (specied as input 
        argument) cloud fraction > 0.7 and cloud tops in the upper 
        troposphere (450-180 hPa).
    (2) Meridional frequency of upper tropospheric (450-180 hPa) cloud 
        fractions > 0.7 for each product.

'''

#    (1) Initial version created on 15 March 2020 by Eloise Marais
#        (eloise.marais@le.ac.uk). Reworked by John Roberts June 2020
#        (jfr10@le.ac.uk)

# Import relevant packages:
import glob
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import argparse
import sys
import os
from os import path

# Import hack
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..'))

from uptrop.gamap_colormap import WhGrYlRd
from uptrop.convert_height_to_press import alt2pres

# Turn off warnings:
# np.warnings.filterwarnings('ignore')

# Define constants/maybe parameters?
FILL_VAL = 9.96921e+36

# Define model grid lat and long limits:
MIN_LAT = -90.
MAX_LAT = 90.
MIN_LON = -180.
MAX_LON = 180.

# Define limits of pressure range of interest (in hPa):
PMIN = 180.
PMAX = 450.

class FileMismatchException(Exception):
    """
    Raised when there are an unequal number of cloud and N02 files
    """
    pass


class LatLonException(Exception):
    """
    Raised when there is a LatLon mismatch
    """
    pass


class ShapeMismatchException(Exception):
    """
    Raised when a FRESCO and DLR file are not the same shape
    """
    pass

class NanInDataException(Exception):
    """
    Raised when filtering has failed
    """


class CloudVariableStore:
    """
    Class containing the running results of various datas for this project; also saving to netCDF and plotting.
    """
    def __init__(self, data_shape):
        """Creates an empty CloudVariableStore of shape data_shape

        :param data_shape: Tuple of (presumably) two elements for shape of cloud data
        :type data_shape: List of length 2
        """
        # TODO Find out from Eloise what these stand for for better variable names
        # Variables explained:
        # KNMI is the cloud product provied with the TROPOMI NO2 data.
        # DLR is the standalone TROPOMI cloud product.
        # Data are: cloud fraction (cf), cloud top (ct), cloud base (cb),
        # optical depth (od; DLR only).        
        self.gknmi_cf = np.zeros(data_shape)  
        self.gknmi_ct = np.zeros(data_shape)
        self.gknmi_cb = np.zeros(data_shape)
        self.gknmi_cnt = np.zeros(data_shape)
        self.gdlr_cf = np.zeros(data_shape)
        self.gdlr_ct = np.zeros(data_shape)
        self.gdlr_cb = np.zeros(data_shape)
        self.gdlr_od = np.zeros(data_shape)
        self.gdlr_cnt = np.zeros(data_shape)
        self.gdiff_cf = np.zeros(data_shape)
        self.gdiff_ct = np.zeros(data_shape)
        self.nobs_dlr = 0
        self.nobs_fresco = 0
        self.shape = data_shape
        # Define data arrays to output cloud fraction binned into 0.05
        # increments from 0.7 to 1.0 centred at 0.725, 0.775 etc.:
        # Consider also splitting this into latitude bands.
        self.cldbin = np.arange(0.7, 1.2 + 0.05, 0.05)
        self.latbin = np.arange(-90, 90 + 15, 15)
        self.knmi_cf_freq = np.zeros((len(self.latbin), len(self.cldbin)))
        self.dlr_cf_freq = np.zeros((len(self.latbin), len(self.cldbin)))

    def update_pixel(self, tropomi_data, trop_i, trop_j):
        """Updates the appropriate pixel of running cloud variables with the tropomi data at trop_i, trop_j.

        :param tropomi_data: The tropomi dataset
        :type tropomi_data: CloudComparisonData
        :param trop_i: First index of pixel
        :param trop_j: Second index of pixel
        :type trop_i: int
        :type trop_j: int
        """
        # TODO: Make this vectorised
        # Skip where FRESCO cloud product is NAN:
        if np.isnan(tropomi_data.tffrc[trop_i, trop_j]) or \
           np.isnan(tropomi_data.tdfrc[trop_i, trop_j]):
            #print("No FRESCO or DLR data, skipping....")
            return

        # Find corresponding gridsquare
        # TODO: Find a better place to put out_lon
        p = np.argmin(abs(out_lon - tropomi_data.tdlons[trop_i, trop_j]))
        q = np.argmin(abs(out_lat - tropomi_data.tdlats[trop_i, trop_j]))

        # Error checks:
        # Check that lat and lon values are the same (except for
        # those at 180/-180 that can have different signs):
        if (tropomi_data.tdlons[trop_i, trop_j] != tropomi_data.tflons[trop_i, trop_j]) \
                and (abs(tropomi_data.tdlons[trop_i, trop_j]) != 180.0):
            print('Longitudes not the same')
            print(tropomi_data.tdlons[trop_i, trop_j], tropomi_data.tflons[trop_i, trop_j])
            print(p)
            print(out_lon[p])
            return
        if tropomi_data.tdlats[trop_i, trop_j] != tropomi_data.tflats[trop_i, trop_j]:
            print('Latitudes not the same')
            return

        # print("{},{} is a valid pixel".format(trop_i, trop_j))

        # Add data to output arrays:
        self.gknmi_cf[p, q] += tropomi_data.tffrc[trop_i, trop_j]
        self.gknmi_ct[p, q] += np.divide(tropomi_data.tftop[trop_i, trop_j], 100)
        self.gknmi_cnt[p, q] += 1.0
        self.gdlr_cf[p, q] += tropomi_data.tdfrc[trop_i, trop_j]
        self.gdlr_ct[p, q] += np.divide(tropomi_data.tdtop[trop_i, trop_j], 100)
        self.gdlr_cb[p, q] += np.divide(tropomi_data.tdbase[trop_i, trop_j], 100)
        self.gdlr_od[p, q] += tropomi_data.tdoptd[trop_i, trop_j]
        self.gdlr_cnt[p, q] += 1.0

    def bin_cloud_fraction(self, tropomi_data):
        """
        # Gather data on frequency of cloud fraction >= 0.7:
        # This is done before filtering for scenes with knmi cloud frac
        # >= 0.7 to also include all relevant DLR scenes:

        :param tropomi_data: The tropomi data containing the cloud fractions
        :type tropomi_data: CloudComparisonData
        """
        # loop over cloud and latitude band bins:
        for w in range(len(self.cldbin)):
            for n in range(len(self.latbin)):

                # Define minimum value in range of interest to avoid cloud fractions < 0.7:
                if w == 0:
                    min_val = self.cldbin[w]
                else:
                    min_val = self.cldbin[w] - 0.025
                

                # (1) KNMI:
                # Get indices for relevant data:
                fbin = np.where((tropomi_data.tffrc >= min_val) &
                                (tropomi_data.tffrc < (self.cldbin[w] + 0.025)) &
                                (tropomi_data.tflats >= (self.latbin[n] - 7.5)) &
                                (tropomi_data.tflats < (self.latbin[n] + 7.5)) &
                                (tropomi_data.tftop >= PMIN*1e2) &
                                (tropomi_data.tftop <= PMAX*1e2) &
                                (~np.isnan(tropomi_data.tffrc)))[0]
                if len(fbin) > 0:
                    self.knmi_cf_freq[n, w] += len(fbin)

                # (2) DLR-OCRA:
                # Get indices for relevant data:
                dbin = np.where((tropomi_data.tdfrc >= min_val) &
                                (tropomi_data.tdfrc < (self.cldbin[w] + 0.025)) &
                                (tropomi_data.tdlats >= (self.latbin[n] - 7.5)) &
                                (tropomi_data.tdlats < (self.latbin[n] + 7.5)) &
                                (tropomi_data.tdtop >= PMIN*1e2) &
                                (tropomi_data.tdtop <= PMAX*1e2) &
                                (~np.isnan(tropomi_data.tdfrc)))[0]
                if len(dbin) > 0:
                    self.dlr_cf_freq[n, w] += len(dbin)

    def update_nobs(self, tropomi_data):
        """Given a tropomi_data object, updates the number of observations.

        :param tropomi_data: The tropomi data objects
        :type tropomi_data: CloudComparisonData
        """
        self.nobs_dlr += tropomi_data.nobs_dlr
        self.nobs_fresco += tropomi_data.nobs_fresco

    def calc_cloud_statistics(self):
        """Finalises cloud statistics

        For each pixel in the final model, calculate the cloud fraction, cloud top pressure, cloud albedo and (for
        DLR) cloud base pressure. Cloud pressure is given in Pa.
        """
        # Get means and differences (only means for cloud base height):
        # (1) Cloud fraction:
        self.gknmi_cf = np.where(self.gknmi_cnt == 0, np.nan, np.divide(self.gknmi_cf, self.gknmi_cnt))
        self.gdlr_cf = np.where(self.gdlr_cnt == 0, np.nan, np.divide(self.gdlr_cf, self.gdlr_cnt))
        self.gdiff_cf = np.subtract(self.gknmi_cf, self.gdlr_cf)
        # (2) Cloud top pressure:
        self.gknmi_ct = np.where(self.gknmi_cnt == 0, np.nan, np.divide(self.gknmi_ct, self.gknmi_cnt))
        self.gdlr_ct = np.where(self.gdlr_cnt == 0, np.nan, np.divide(self.gdlr_ct, self.gdlr_cnt))
        self.gdiff_ct = np.subtract(self.gknmi_ct, self.gdlr_ct)
        # (3) Cloud albedo/optical depth:
        self.gdlr_od = np.where(self.gdlr_cnt == 0, np.nan, np.divide(self.gdlr_od, self.gdlr_cnt))
        # (4) Cloud base pressure (DLR only):
        self.gdlr_cb = np.divide(self.gdlr_cb, self.gdlr_cnt)
        self.gdlr_cb[self.gdlr_cnt == 0.0] = np.nan

    def write_to_netcdf(self, out_dir):
        """Given an out_directory, writes totalled data to netcdf.

        The file will be named 'fresco-dlr-cloud-products-[MMName]-[StrYY]-[out_res]-ref-[ref_cld_prod]-v3.nc'

        :param out_dir: The directory that will contain the file.
        :type out_dir: str
        """
        out_dir = path.abspath(out_dir)

        # Write data to file:
        outfile = path.join(out_dir, 'fresco-dlr-cloud-products-' + MMName + '-' + StrYY + '-' + out_res + '-ref-' + ref_cld_prod + '-v3.nc')
        ncfile = Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
        ncfile.createDimension('xdim', len(out_lon))    # TODO: Make this a member
        ncfile.createDimension('ydim', len(out_lat))
        ncfile.createDimension('frdim', len(self.cldbin))
        ncfile.createDimension('lbdim', len(self.latbin))
        # Global attributes:
        ncfile.title = 'FRESCO and DLR TROPOMI monthly mean cloud properties for ' + \
                       StrMM + ' ' + StrYY
        ncfile.subtitle = 'Data written to file in ' + StrMM + ' ' + StrYY
        ncfile.anything = 'Verions used are v010107 for the DLR product and v010302 for the FRESCO product'
        # Longitudes:
        lon = ncfile.createVariable('lon', np.dtype('float'), ('xdim',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        lon[:] = out_lon
        # Latitudes:
        lat = ncfile.createVariable('lat', np.dtype('float'), ('ydim',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lat[:] = out_lat
        # Cloud fraction bins:
        cbin = ncfile.createVariable('cld_fr_bin', np.dtype('float'), ('frdim',))
        cbin.units = 'unitless'
        cbin.long_name = 'Cloud fraction bin centre'
        cbin[:] = self.cldbin

        # Latitude bands bins:
        lbbin = ncfile.createVariable('lat_band_bin', np.dtype('float'), ('lbdim',))
        lbbin.units = 'unitless'
        lbbin.long_name = 'Latitude band bin centre'
        lbbin[:] = self.latbin
        # DLR cloud fraction:
        data1 = ncfile.createVariable('dlr_cld_frac', np.dtype('float'), ('xdim', 'ydim'))
        data1.units = 'unitless'
        data1.long_name = 'DLR CAL cloud fraction'
        data1[:] = self.gdlr_cf
        # DLR cloud top pressure:
        data2 = ncfile.createVariable('dlr_cld_top_pres', np.dtype('float'), ('xdim', 'ydim'))
        data2.units = 'hPa'
        data2.long_name = 'DLR CAL cloud top pressure'
        data2[:] = self.gdlr_ct
        # DLR cloud base pressure:
        data3 = ncfile.createVariable('dlr_cld_base_pres', np.dtype('float'), ('xdim', 'ydim'))
        data3.units = 'hPa'
        data3.long_name = 'DLR CAL cloud base pressure'
        data3[:] = self.gdlr_cb
        # DLR cloud optical depth:
        data4 = ncfile.createVariable('dlr_cld_optd', np.dtype('float'), ('xdim', 'ydim'))
        data4.units = 'unitless'
        data4.long_name = 'DLR CAL cloud optical thickness'
        data4[:] = self.gdlr_od
        # FRESCO cloud fraction:
        data5 = ncfile.createVariable('knmi_cld_frac', np.dtype('float'), ('xdim', 'ydim'))
        data5.units = 'unitless'
        data5.long_name = 'FRESCO cloud fraction'
        data5[:] = self.gknmi_cf
        # FRESCO cloud top pressure:
        data6 = ncfile.createVariable('knmi_cld_top_pres', np.dtype('float'), ('xdim', 'ydim'))
        data6.units = 'hPa'
        data6.long_name = 'FRESCO cloud top pressure'
        data6[:] = self.gknmi_ct
        # FRESCO data points in each gridsquare:
        data7 = ncfile.createVariable('knmi_points', np.dtype('float'), ('xdim', 'ydim'))
        data7.units = 'unitless'
        data7.long_name = 'Number of FRESCO data points'
        data7[:] = self.gknmi_cnt
        # DLR data points in each gridsquare:
        data8 = ncfile.createVariable('dlr_points', np.dtype('float'), ('xdim', 'ydim'))
        data8.units = 'unitless'
        data8.long_name = 'Number of DLR data points'
        data8[:] = self.gdlr_cnt
        # Frequency of FRESCO cloud fractions > 0.7:
        data9 = ncfile.createVariable('knmi_cf_freq', np.dtype('float'), \
                                      ('lbdim', 'frdim',))
        data9.units = 'unitless'
        data9.long_name = 'Freqeuncy distribution of FRESCO cloud fractions > 0.7'
        data9[:] = self.knmi_cf_freq
        # Frequency of DRL-OCRA cloud fractions > 0.7:
        data10 = ncfile.createVariable('dlr_cf_freq', np.dtype('float'), ('lbdim', 'frdim'))
        data10.units = 'unitless'
        data10.long_name = 'Freqeuncy distribution of DLR cloud fractions > 0.7'
        data10[:] = self.dlr_cf_freq

        # close the file.
        ncfile.close()

    def plot_clouds_products(self, plot_dir):
        """For cloud fraction, cloud top pressure, plots the products for DLR, FRESCO and DLR-FRESCO

        :param plot_dir: The directory that will contain the plots
        :type plot_dir: str"""
        # TODO: Get MMName and StrYY from file_path instead of the global at the bottom of the script
        # PLOT THE DATA:
        m = Basemap(resolution='l', projection='merc', lat_0=0, lon_0=0,
                    llcrnrlon=MIN_LON, llcrnrlat=-70,
                    urcrnrlon=MAX_LON, urcrnrlat=70)
        xi, yi = m(X, Y)
        # (1) Cloud fraction:
        plt.figure(1)
        plt.subplot(3, 1, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.gknmi_cf), vmin=0.7, vmax=1.3, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('FRESCO Cloud fraction')
        plt.subplot(3, 1, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_cf), vmin=0.7, vmax=1.3, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR Cloud fraction')
        plt.subplot(3, 1, 3)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdiff_cf), vmin=-0.2, vmax=0.2, cmap='bwr')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('Difference (FRESCO-DLR)')
        plt.savefig(path.join(plot_dir, 'fresco-vs-dlr-cloud-frac-' + MMName + '-' + StrYY + '.ps'), \
                    format='ps')
        # (2) Cloud top pressure:
        plt.figure(2)
        plt.subplot(3, 1, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.gknmi_ct), vmin=150, vmax=500, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('FRESCO Cloud top pressure [hPa]')
        plt.subplot(3, 1, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_ct), vmin=150, vmax=500, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR Cloud top pressure [hPa]')
        plt.subplot(3, 1, 3)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdiff_ct), vmin=-30, vmax=30, cmap='bwr')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('Difference (FRESCO-DLR)')
        plt.savefig(path.join(plot_dir, 'fresco-vs-dlr-cloud-top-press-' + MMName + '-' + StrYY + '.ps'),
                    format='ps')
        # (3) Cloud optical depth/albedo (DLR only):
        plt.figure(3)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_od), vmin=0, vmax=200, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR cloud optical thickness')
        plt.savefig(path.join(plot_dir, 'dlr-cloud-optical-depth-' + MMName + '-' + StrYY + '.ps'),
                    format='ps')
        # (4) Cloud base pressure (DLR only):
        plt.figure(4)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_cb), vmin=150, vmax=500, cmap='jet')
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR base pressure [hPa]')
        plt.savefig(path.join(plot_dir,'dlr-cloud-base-press-' + MMName + '-' + StrYY + '.ps'),
                    format='ps')
        # (5) Number of points (both):
        plt.figure(4)
        plt.subplot(2, 1, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.gknmi_cnt), vmin=0, vmax=500, cmap=WhGrYlRd)
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('FRESCO No. of obs')
        plt.subplot(2, 1, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.gdlr_cnt), vmin=0, vmax=500, cmap=WhGrYlRd)
        m.drawcoastlines()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('DLR No. of obs')
        plt.savefig(path.join(plot_dir, 'fresco-dlr-number-of-obs-' + MMName + '-' + StrYY + '.ps'),
                    format='ps')
        plt.show()


class CloudComparisonData:
    """
    Class for holding the data for a fresco-dlr file pair. Applies data filters on creation.
    """

    def __init__(self, td_file_path, tf_file_path, ref_cld_prd):
        """Returns an instance of CloudComparisonData for comparing fresco and dlr files

        :param td_file_path: The path to the DLR file
        :type td_file_path: str
        :param tf_file_path: The path to the Fresco file
        :type tf_file_path: str
        :param ref_cld_prd: The reference cloud product for identifying optically thick clouds in the upper troposphere
        :type ref_cld_prd: str
        :return: A filtered and sanity-checked cloud comparison
        :rtype: CloudComparisonData
        """
        # Extract orbit numbers from the file name.
        self.forb = path.basename(tf_file_path)[104:109]
        self.dorb = path.basename(td_file_path)[106:111]
        # Check orbit/swath is the same. If not, skip this iteration:
        if (self.forb != self.dorb):
            print("Orbit is not swath")
            return
        self.read_tdfile(td_file_path)
        self.read_tffile(tf_file_path)
        self.ref_cld_prd = ref_cld_prd
        self.nobs_dlr, self.nobs_fresco = self.get_nobs()
        self.filter_tdfile()
        self.filter_tffile()
        self.check_parity()
        self.shape = self.tdlons.shape

    def read_tdfile(self, tdfile):
        """Read DLR cloud data

        :param tdfile: Path to the DLR file
        :type tdfile: str"""
        dlr_cloud_data = Dataset(tdfile, mode='r')
        # Extract data of interest:
        self.tdlons = dlr_cloud_data.groups['PRODUCT'].variables['longitude'][:].data[0, :, :]
        self.tdlats = dlr_cloud_data.groups['PRODUCT'].variables['latitude'][:].data[0, :, :]
        tfrc = dlr_cloud_data.groups['PRODUCT'].variables['cloud_fraction'][:]
        self.tdfrc = tfrc.data[0, :, :]
        if ( dlr_cld_top=='height' ):
            ttop = dlr_cloud_data.groups['PRODUCT'].variables['cloud_top_height'][:]
            self.tdtop = ttop.data[0, :, :]       
            # Convert height to pressure (done before counting the number
            # if valid data points in get_nobs:            
            self.tdtop = np.where(self.tdtop != FILL_VAL, alt2pres(self.tdtop), FILL_VAL) 
        else:
            ttop = dlr_cloud_data.groups['PRODUCT'].variables['cloud_top_pressure'][:]     
            self.tdtop = ttop.data[0, :, :]
        tbase = dlr_cloud_data.groups['PRODUCT'].variables['cloud_base_pressure'][:]
        self.tdbase = tbase.data[0, :, :]
        tqval = dlr_cloud_data.groups['PRODUCT'].variables['qa_value'][:]
        self.tdqval = tqval.data[0, :, :]  # for more accurate CAL product
        toptd = dlr_cloud_data.groups['PRODUCT'].variables['cloud_optical_thickness'][:]
        self.tdoptd = toptd.data[0, :, :]
        tsnow = dlr_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['snow_ice_flag'][:]
        self.tdsnow = tsnow.data[0, :, :]
        dlr_cloud_data.close()   # Note to self; does NetCDF have a context manager?

    def read_tffile(self, tffile):
        """Read in FRESCO cloud data

        :param tffile: Path to the tffile
        :type tffile: str"""
        fresco_cloud_data = Dataset(tffile, mode='r')
        # Extract data of interest:
        self.tflons = fresco_cloud_data.groups['PRODUCT'].variables['longitude'][:].data[0, :, :]
        self.tflats = fresco_cloud_data.groups['PRODUCT'].variables['latitude'][:].data[0, :, :]
        tfrc = fresco_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['cloud_fraction_crb'][:]
        self.tffrc = tfrc.data[0, :, :]
        talb = fresco_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['cloud_albedo_crb'][:]
        self.tfalb = talb.data[0, :, :]
        ttop = fresco_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['cloud_pressure_crb'][:]
        self.tftop = ttop.data[0, :, :]
        tqval = fresco_cloud_data.groups['PRODUCT'].variables['qa_value'][:]
        self.tfqval = tqval.data[0, :, :]
        tsnow = fresco_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['snow_ice_flag'][:]
        self.tfsnow = tsnow.data[0, :, :]
        tscenep = fresco_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['apparent_scene_pressure'][:]
        self.tfscenep = tscenep.data[0, :, :]
        tsurfp = fresco_cloud_data.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['surface_pressure'][:]
        self.tfsurfp = tsurfp.data[0, :, :]
        fresco_cloud_data.close()

    def filter_tdfile(self):
        """Filters DLR data 

        Applies the following filters to the dlr data
        Snow: Drops values flagged as snow
        Poor quality data: quality flag is < 0.5
        Drop the following values if DLR-OCRA is the reference cloud product:
         - cloud fraction < 0.7
         - cloud top pressure is between PMAX and PMIN

        :raises NanInDataException: Raised when any NaNs remain in data
        """

        # Set missing/poor quality/irrelevant data to NAN:
        # Apply cloud fraction filter, but set it to 0.7 rather
        #     than 0.9 to account for variability around this threshold
        #     in the two cloud products.
        if self.ref_cld_prd=='dlr-ocra':
            self.tdfrc = np.where(self.tdfrc < 0.7, np.nan, self.tdfrc)
            self.tdtop = np.where(self.tdtop > (PMAX * 1e2), np.nan, self.tdtop)
            self.tdtop = np.where(self.tdtop < (PMIN * 1e2), np.nan, self.tdtop)
        # Fill values (do for cloud fraction and cloud top pressure, as missing
        # values for cloud fraction may not be the same as missing values for
        # other data, as these are obtained with different algorithms):
        # necessarily the same as the fill values for the other data):
        self.tdfrc = np.where(self.tdfrc == FILL_VAL, np.nan, self.tdfrc)
        self.tdtop = np.where(self.tdtop == FILL_VAL, np.nan, self.tdtop)
        # Set cloud fraction to nan for scenes with missing cloud top
        # pressure data:
        self.tdfrc = np.where(np.isnan(self.tdtop), np.nan, self.tdfrc)
        # Poor data quality:
        self.tdfrc = np.where(self.tdqval < 0.5, np.nan, self.tdfrc)
        # Snow/ice cover:
        self.tdfrc = np.where(self.tdsnow != 0, np.nan, self.tdfrc)

        # Apply filter to remaining data:
        self.tdtop = np.where(np.isnan(self.tdfrc), np.nan, self.tdtop)
        self.tdbase = np.where(np.isnan(self.tdfrc), np.nan, self.tdbase)
        self.tdoptd = np.where(np.isnan(self.tdfrc), np.nan, self.tdoptd)
        # Error check:
        if np.nanmax(self.tdtop) == FILL_VAL:
            raise NanInDataException('Not all missing values converted to NANs')
        if np.nanmax(self.tdbase) == FILL_VAL:
            raise NanInDataException('Not all missing values converted to NANs')
        if np.nanmax(self.tdoptd) == FILL_VAL:
            raise NanInDataException('Not all missing values converted to NANs')

    def filter_tffile(self):
        """Filters fresco data

        Applies the following filters to the fresco data
        Snow: Drops values for snow, coastline, snow/ice misclassified as clouds
        Poor quality data: quality flag is < 0.45
        Drop the following values if FRESCO is the reference cloud product:
         - cloud fraction < 0.7
         - cloud top pressure is between PMAX and PMIN
        """
        # Convert all valid snow/ice free flag values (0,255) to 0.
        self.tfsnow = np.where(self.tfsnow == 255, 0, self.tfsnow)
        # Coastlines (listed as potential "suspect" in the ATBD document p. 67):
        self.tfsnow = np.where(self.tfsnow == 252, 0, self.tfsnow)
        # Less then 1% snow/ice cover:
        self.tfsnow = np.where(self.tfsnow < 1, 0, self.tfsnow)
        # Snow/ice misclassified as clouds:
        self.tfsnow = np.where((self.tfsnow > 80)
                               & (self.tfsnow < 104)
                               & (self.tfscenep > (np.multiply(0.98, self.tfsurfp))),
                               0, self.tfsnow)

        # Set missing/poor quality/irrelevant data to NAN:
        # Apply cloud fraction filter, but set it to 0.7 rather
        #     than 0.9 to account for variability around this threshold
        #     in the two cloud products.
        if self.ref_cld_prd=='fresco':
            self.tffrc = np.where(self.tffrc < 0.7, np.nan, self.tffrc)
            # Apply cloud top pressure filter:
            self.tffrc = np.where(self.tftop > (PMAX * 1e2), np.nan, self.tffrc)
            self.tffrc = np.where(self.tftop < (PMIN * 1e2), np.nan, self.tffrc)
        self.tffrc = np.where(self.tffrc == FILL_VAL, np.nan, self.tffrc)
        # QA Flags. Threshold of 0.45 suggested by Henk Eskes in email
        # exchange on 18 Jan 2020:
        self.tffrc = np.where(self.tfqval < 0.45, np.nan, self.tffrc)
        # Snow/ice cover:
        self.tffrc = np.where(self.tfsnow != 0, np.nan, self.tffrc)
        # Apply filter to remaining data:
        self.tftop = np.where(np.isnan(self.tffrc), np.nan, self.tftop)

    def check_parity(self):
        """Checks the shape of fresco and dlr data match

        :raises ShapeMismatchException: Raised when the shapes do not match"""
        # Skip files if the number of indices are not equal:
        if self.tdlons.shape != self.tflons.shape:
            print('Indices not equal for orbit {}'.format(self.forb))
            raise ShapeMismatchException

    def get_nobs(self):
        """Gets the number of valid observations in this cloud data

        Returns a tuple of number of valid cloud observations of (dlr product, fresco product)
        A valid observation has:
         - Quality of <0.5 (dlr) or <0.45(fresco)
         - Cloud fraction >= 0.7
         - top-of-atmosphere pressure between PMIN and PMAX
         - No snow (snow mask = 0)

        :return: A tuple of number of valid observations of (dlr, fresco)
        :rtype: tuple (int, int)
        """
        # Bug might still be here? I think we fixed it.
        dlr_ind = np.count_nonzero((self.tdqval < 0.5) & (self.tdfrc >= 0.7) & (self.tdtop >= PMIN*1e2) & (self.tdtop <= PMAX*1e2) & (self.tdsnow == 0))
        fr_ind = np.count_nonzero((self.tfqval < 0.45) & (self.tffrc >= 0.7) & (self.tftop >= PMIN*1e2) & (self.tftop <= PMAX*1e2) & (self.tfsnow == 0))
        # DLR value will be different, fr_ind will be same
        return dlr_ind, fr_ind


def process_file(tdfile, tffile, running_total_container):
    """Processes a paired dlr and fresco product, adds the pixels to the running total in the
    running_total_container and updates the cloud_fraction

    :param tdfile: Path to the dlr file
    :type tdfile: str
    :param tffile: Path to the fresco file
    :type tffile: str
    :param running_total_container: Instance of CloudVariableStore to update with the dlr and fresco file
    :type running_total_container: CloudVariableStore
    """
    # Track progress:
    print('===> Processing: ', tdfile)
    try:
        file_data_container = CloudComparisonData(tdfile, tffile, ref_cld_prod)
        running_total_container.update_nobs(file_data_container)
        #print("Fresco nobs: {}\nDLR nobs: {}".format(
        #    file_data_container.nobs_fresco, file_data_container.nobs_dlr))
        running_total_container.bin_cloud_fraction(file_data_container)
        # REGRID THE DATA:
        for i in range(file_data_container.shape[0]):
            for j in range(file_data_container.shape[1]):
                running_total_container.update_pixel(file_data_container, i, j)
    except ShapeMismatchException:
        print("Mismatch in shape of {} and {}".format(tdfile, tffile))
    return


def get_files_for_month(sen_5_p_dir, month_index, ndays=31):
    """Gets fresco an
    For a given month index (jan-may 2020 being 1-5, jun-dec 2019 being 6-12), returns every DLR and Fresco
    file_path. Also sets the globals StrMM, StrYY, MMName (used in the plotting method)(at least until I fix it)

    :param sen_5_p_dir: The directory containing the DLR and Fresco files
    :type sen_5_p_dir: str
    :param month_index: The month index
    :type month_index: int
    :param ndays: The number of days in the month
    :type ndays: int

    :raises FileMismatchException: Raised if there is an unequal number of fresco and dlr files

    :return: The list of dlr files and the list of fresco files
    :rtype: tuple(dlr files, fresco files)
    """
    # TODO Roll the string manufacturing into the CloudVariableStore class
    global StrMM, StrYY, MMName


    # Input parameter (to selet month of interest):
    # Define month of interest as string 2 characters in length:
    # TODO: Change this entire section into datetime
    StrMM = str(month_index)
    # Define string of year and first 3 letters of month name based on above entry:
    if StrMM == '01': StrYY, MMName = '2020', 'jan'
    if StrMM == '02': StrYY, MMName = '2020', 'feb'
    if StrMM == '03': StrYY, MMName = '2020', 'mar'
    if StrMM == '04': StrYY, MMName = '2020', 'apr'
    if StrMM == '05': StrYY, MMName = '2020', 'may'
    if StrMM == '06': StrYY, MMName = '2019', 'jun'
    if StrMM == '07': StrYY, MMName = '2019', 'jul'
    if StrMM == '08': StrYY, MMName = '2019', 'aug'
    if StrMM == '09': StrYY, MMName = '2019', 'sep'
    if StrMM == '10': StrYY, MMName = '2019', 'oct'
    if StrMM == '11': StrYY, MMName = '2019', 'nov'
    if StrMM == '12': StrYY, MMName = '2019', 'dec'
    # Define days:
    dd = list(range(1, ndays + 1))
    strdd = [''] * ndays
    cnt = 0
    for d in dd:
        strdd[cnt] = '0' + str(d) if d < 10 else str(d)
        cnt = cnt + 1

    # dir structure specific to HPC, file_path from TROPOMI
    # Get DLR data file names:
    td_file_list = glob.glob(
        path.join(sen_5_p_dir,
                  'CLOUD_OFFL',
                  StrYY,
                  StrMM,
                  'S5P_OFFL_L2__CLOUD__' + StrYY + StrMM + '*'))
    td_file_list = sorted(td_file_list)
    # Get FRESCO file names:
    tf_file_list = glob.glob(
        path.join(sen_5_p_dir,
                  'NO2_OFFL',
                  StrYY,
                  StrMM ,
                  'S5P_OFFL_L2__NO2____' + StrYY + StrMM + '*'))
    tf_file_list = sorted(tf_file_list)
    # Check that number of files are equal. If not, exit the programme:
    if len(td_file_list) != len(tf_file_list):
        print('DLR files = ', len(td_file_list))
        print('FRESCO files = ', len(tf_file_list))
        print('unequal number of files')
        raise FileMismatchException("Unequal number of DLR and FRESCO files, check archive.")
    if len(td_file_list) == 0:
        raise FileNotFoundError("No files found in {}".format(sen_5_p_dir))

    return td_file_list, tf_file_list


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Extracts and plots cloud data")
    parser.add_argument("--s5p_data_dir", default='/data/uptrop/nobackup/tropomi/Data/')
    parser.add_argument("--output_dir", default='/data/uptrop/Projects/UpTrop/python/Data_v2/')
    parser.add_argument("--month", default="10")
    parser.add_argument("--plot_dir", default="/data/uptrop/Projects/UpTrop/python/Images/")
    parser.add_argument("--number_of_days", default=31, type=int)
    parser.add_argument("--out_res", default="1x1")
    parser.add_argument("--dlr_cld_top", default="height")
    parser.add_argument("--ref_cld_prod", default="fresco")
    args = parser.parse_args()

    s5p_data_dir = path.expanduser(args.s5p_data_dir)
    output_dir = path.expanduser(args.output_dir)
    plot_dir = path.expanduser(args.plot_dir)
    dlr_cld_top = path.expanduser(args.dlr_cld_top)
    ref_cld_prod = path.expanduser(args.ref_cld_prod)

    # TODO: Make this a member of CloudVariableStore
    if args.out_res == '1x1':
        DELTA_LAT = 1
        DELTA_LON = 1
    if args.out_res == '2x25':
        DELTA_LAT = 2
        DELTA_LON = 2.5
    if args.out_res == '4x5':
        DELTA_LAT = 4
        DELTA_LON = 5
    out_lon = np.arange(MIN_LON, MAX_LON+DELTA_LON, DELTA_LON)
    out_lat = np.arange(MIN_LAT, MAX_LAT+DELTA_LAT, DELTA_LAT)
    # Convert output lats and long to 2D:
    X, Y = np.meshgrid(out_lon, out_lat, indexing='ij')

    out_res = path.expanduser(args.out_res)

    td_file_list, tf_file_list = get_files_for_month(s5p_data_dir, args.month, args.number_of_days)
    running_cloud_total = CloudVariableStore(X.shape)

    # Loop over files:
    for td_file, tf_file in zip(td_file_list, tf_file_list):
        print("Processing files {} and {}".format(tf_file, td_file))
        process_file(td_file, tf_file, running_cloud_total)

    running_cloud_total.calc_cloud_statistics()

    # Print number of observations to screen:
    print('No. of FRESCO obs for '+MMName+' = ',running_cloud_total.nobs_fresco)
    print('No. of DLR obs for '+MMName+' = ',running_cloud_total.nobs_dlr)
    print("Writing to NetCDF at {}".format(args.output_dir))
    running_cloud_total.write_to_netcdf(output_dir)
    print("Creating plots at {}".format(args.plot_dir))
    running_cloud_total.plot_clouds_products(plot_dir)
    print("Processing complete.")

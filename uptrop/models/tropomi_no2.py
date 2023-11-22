import os
import logging
from netCDF4 import Dataset
import numpy as np
from ..utils.dates_files_utils import get_date


class TropomiNO2:
    """A class for extracting, preprocessing and containing data from a s5p tropomi file."""

    # Initialize _everything_ before the first read.
    def __init__(self, file_path, no2_prod):
        """Reads the tropomi file at file_path and prepares it for processing.

        :param file_path: Path to the netcdf4 file containing the tropomi data
        :type file_path: str
        """

        self.file_name = os.path.basename(file_path)
        logging.info(f'Processing: {self.file_name}')
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
            error_message = (
                f'NO2 file: {self.date}, Cloud file: {cloud_data.date}.\n'
                f'EXITING: Files are not for the same date!'
            )
            logging.info(error_message)
            raise ValueError(error_message)
    
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
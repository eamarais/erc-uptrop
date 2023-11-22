import os
import logging
import numpy as np
from netCDF4 import Dataset
from .helpers import alt2pres
from ..utils.dates_files_utils import get_date


class CloudDataTropomiNO2:
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
        self.file_name = os.path.basename(file_path)
        self.date = get_date(self.file_name)  # Assuming for now that ocra and S5P have the same timestamping

        # Set from input:
        self.fresco_440 = fresco_440

        # Set from file contents
        self.cldfrac = None
        self.tcldpres = None
        self.tsnow = None

        # Get cloud fields: 

        """ 
        Need to add a check here, if the combo of cloud and TROPOMI products are not available,
        terminal the programme and notify the user.
        """
        valid_combinations = [
        ('fresco-wide', 'PAL'),
        ('fresco-wide', 'OFFL'),
        ('o22cld', 'PAL'),
        ('dlr-ocra', 'OFFL'),
        ]

        if (data_type, no2_prod) not in valid_combinations:
            print("Invalid combination of data_type and no2_prod. Program terminated.")
            return
    
        if data_type == 'fresco-wide' and no2_prod == 'PAL':
            self.get_fresco_cloud_fields(file_path)
            # Fix snow/ice/coast flag issue in PAL_ files:
            self.fix_snow()
        elif data_type == 'fresco-wide' and no2_prod == 'OFFL':
            self.get_fresco_cloud_fields(file_path)
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
            logging.info('this method of defining the fill value for dlr-ocra does not work. FIX!!!')

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
        #    logging.info('this method of defining the fill value for dlr-ocra does not work. FIX!!!')
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

    def check_parity(self, trop_data):
        # Skip files if the number of indices are not equal:            
        if self.cldfrac.shape != trop_data.sza.shape:
            logging.info(f'Cloud product and NO2 indices ne!')
            logging.info(f'{self.cldfrac.shape}, {trop_data.sza.shape}')
            logging.info(f'Skipping this swath')
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
        
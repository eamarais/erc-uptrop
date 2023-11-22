# GridAggregator for Tropomi O3

import logging
import numpy as np
from netCDF4 import Dataset
from .helpers import cloud_slice, CLOUD_SLICE_ERROR_ENUM


class UnequalColumnException(Exception):
    pass


class GridAggregatorTropomiO3:
    """A class for aggregating higher-resolution data into grid squares"""
    def __init__(self, dellat, dellon, pmin, pmax, first):
        """Creates a grid aggregator across the entire world of resolution dellon, dellat

        :param dellat: vertical size of the aggregation grid in decimal degrees
        :type dellat: float
        :param dellon: Horizontal size of the aggregation grid in decimal degrees
        :type dellon: float

        :returns: A GridAggregator object
        :rtype: GridAggregator"""
        if first==0:
            self.postfilt = []
            self.out_lon = np.arange(-180, 180 + dellon, dellon)
            self.out_lat = np.arange(-90, 90 + dellat, dellat)

            self.xdim = len(self.out_lon)
            self.ydim = len(self.out_lat)

            # Seasonal mean output arrays
            self.go3vmr = np.zeros((self.xdim, self.ydim))  # o3 VMR in pptv
            self.gcnt = np.zeros((self.xdim, self.ydim))  # No. of data points (orbits)
            self.gerr = np.zeros((self.xdim, self.ydim))  # Weighted error
            self.gwgt = np.zeros((self.xdim, self.ydim))  # Weights
            self.pcld_range = np.zeros((self.xdim, self.ydim))  # Cloud top pressure range used for cloud-slicing
            self.pcld_ceil = np.zeros((self.xdim, self.ydim))  # Ceiling (min pressure) of the cloud top pressure range
            
            self.file_count = 0
            self.current_max_points = 0

            # Calculate pressure range and centre for Gaussian weighting:
            self.pmin = pmin
            self.pmax = pmax
            self.press_sigma = 0.5*(pmax - pmin)
            self.press_mid   = pmin + self.press_sigma

            # Define cloud pressure difference threshold:
            # Use half the range cloud height range so that code is versatile and
            # because this is consistent with Marais et al. (2021) for the upper troposphere.
            self.diff_cldh_thold = self.press_sigma*1.2
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
                "o3_outlier": 0,
                "non_uni_strat": 0,
            }

            # Members from add_trop_data_to_gridsquare
            self.go3 = None
            self.gstrat = None
            self.gcldp = None
            self.cntloop = None

            self.cloud_slice_count = 0
            first = 1

    def initialise_grid(self):
        """Zeros the go3, gstrat, gcldp and cntloop members with empty appendable lists"""
        self.go3 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.gstrat = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.gcldp = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.cntloop = [[0 for n in range(self.ydim)] for m in range(self.xdim)]

    def grid_trop_data(self, trop_data):
        """Allocates the strat, gc_data and cloud pressure in trop_data into the go3, gstrat and gcldp grid

        :param trop_data: an instance of TropomiData
        :type trop_data: uptrop.tropomi_ut.gc_data.TropomiData
        """

        for tot_value, trop_lat, trop_lon, cloud_pressure \
         in zip(trop_data.totvcd, trop_data.lats, trop_data.lons, trop_data.cldpres):
            # Skip over pixels where total column is less than stratospheric
            # column. This addresses positive bias in the cloud-sliced results
            # at low concentrations of UT o3:
            #if geo_tot_value < strat_o3_val:
            #    continue
            
            # Find nearest gridsquare in output grid:
            p = np.argmin(abs(self.out_lon - trop_lon))
            q = np.argmin(abs(self.out_lat - trop_lat))            

            # Convert columns from mol/m2 to molec/cm2:
            tvcdo3 = np.multiply( tot_value, trop_data.o3sfac ) #np.multiply(geo_tot_value, trop_data.o3sfac)
            #tstrat = np.multiply(strat_o3_val, trop_data.o3sfac)

            # Skip if o3 value is infinite:
            #if math.isinf(tvcdo3)==True:
            #    logging.info(tot_value)
            #    sys.exit()

            # Get relevant data in each output grid square:
            self.go3[p][q].append(tvcdo3)
            #self.gstrat[p][q].append(tstrat)
            self.gcldp[p][q].append(cloud_pressure)

            # Increment indices:
            self.cntloop[p][q] += 1

        # Save % valid TROPOMI pixels retained before applying cloud
        # slicing:
        self.postfilt.append(100. * (trop_data.tcnt / trop_data.inicnt))
        self.file_count += 1

    def _add_slice(self, i, j, t_cld, t_col_o3):
        """Extracts the upper troposphere gc_data, gc_data error and mean cloud pressure for grid square [i,j]

        This method uses the cloud-slicing function :ref:`uptrop.cloud_slice_ut_o3.cldslice`
        Once calculated, the a weighting is derived from cloud pressure.
        The weighted upper tropospheric gc_data and error is added to the rolling total for this season.
        If the cloud slicing fails, then the reason is added to loss_count for the end report.

        :param i: X-index of grid square
        :type i: int
        :param j: Y-index of grid square
        :type j: int
        :param t_cld: A list of cloud pressures
        :type t_cld: list of floats
        :param t_col_o3: A list of gc_data values, of same length as t_cld
        :type t_col_o3: list of floats
        """
        # Check here, if the code works, is your cloud_slice_o3 wrong?
        # Or is                                     
        utmro3, utmro3err, stage_reached, mean_cld_pres = cloud_slice(
            species="o3",
            pcolo3 = t_col_o3,
            cldtophgt = t_cld,
            cld_diff_thold = self.diff_cldh_thold,
            pmin=self.pmin,
            pmax=self.pmax
            )
        # Calculate weights:
        #gaus_wgt = np.exp((-(mean_cld_pres - 315) ** 2) / (2 * 135 ** 2))
        # Skip if approach didn't work (i.e. cloud-sliced UT o3 is NaN):
        # Drop out after the reason for data loss is added to loss_count.
        if np.isnan(utmro3) or np.isnan(utmro3err):
            # Cloud-slicing led to nan, but due to rma regression rather
            # than data filtering (these are rare):
            if (stage_reached == 0):
                logging.info(f"Cloud-sliced o3 NAN for pixel {i} {j}")
                return
            self.loss_count[CLOUD_SLICE_ERROR_ENUM[stage_reached]] += 1
            # logging.info("Cloud-slice exception {} in pixel i:{} j:{}".format(
            #    CLOUD_SLICE_ERROR_ENUM[stage_reached], i, j))
        else:
            # Calculate weights:
            gaus_wgt = np.exp((-(mean_cld_pres - self.press_mid) ** 2) / (2 * self.press_sigma ** 2))
            #gaus_wgt = (1 / np.square(utmro3err) )
            self.go3vmr[i, j] += np.multiply(utmro3, gaus_wgt)
            #self.go3vmr[i, j] += utmro3 #/ np.square(utmro3err)
            self.pcld_range[i, j] += np.nanmax(t_cld) - np.nanmin(t_cld)
            self.pcld_ceil[i, j] += np.nanmin(t_cld)
            self.gwgt[i, j] += gaus_wgt
            self.gerr[i, j] += utmro3err
            self.gcnt[i, j] += 1
            self.cloud_slice_count += 1

    def apply_cloud_slice(self, n_slices=40):
        """Applies the cloud slicing algorithm to each square in the populated grid.

        This method walks over every gridsquare in go3/gstrat/gcldp.
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

                tcolo3 = self.go3[i][j]
                #strat = self.gstrat[i][j]
                tcld = self.gcldp[i][j]

                # Only loop over grids with relevant data (identified as
                # non-empty lists):
                #if not tcolo3:
                #    continue

                # Skip if fewer than 20 points:
                #if len(tcolo3) < 10:
                #    self.loss_count["too_few_points"] += 1
                #    continue

                # Convert from Pa to hPa for input to the cloud-slicing algorithm:
                #tcolo3 = np.multiply(tcolo3, 1e4)
                tcld = np.multiply(tcld, 1e-2)
                #strat = np.array(strat)

                # Error check that the sizes of the arrays are equal:
                if (len(tcld) != len(tcolo3)):
                    logging.info('Array sizes ne: cloud height and partial column')
                    raise UnequalColumnException

                # Skip scenes with non-uniform stratosphere using the
                # same threshold as is used for GEOS-Chem:
                #if (np.std(strat) / np.mean(strat)) > 0.02:
                #    self.loss_count["non_uni_strat"] += 1
                #    continue

                # Filter stratospheric variability for ozone (lines above are for o3).
                # Use total column ozone where most of ozone is located to
                # identify variability in stratospheric ozone following approach in
                # Ziemke et al. (2003), 
                # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2002JD002919
                # where scenes with horizontal variability < 50 DU are retained
                #if np.abs(np.nanmax(tcolo3) - np.nanmin(tcolo3)) >= 10:
                #    self.loss_count["non_uni_strat"] += 1
                #    continue

                # Get number of points:
                npnts = len(tcld)
                if npnts > self.current_max_points:
                    self.current_max_points = npnts
                    logging.info(f'{self.current_max_points}')

                # Use cloud_slice_ut_o3 function to get o3 mixing
                # ratio from cloud-slicing:
                # Change this slightly to exploit higher spatial resolution of TROPOMI to
                # increase the number of scenes retrieved:
                #if ((npnts >= 10) & (npnts < 50)):
                self._add_slice(i,j,tcld,tcolo3)

                #elif (npnts >= 50):
                #    # Define number of iterations:
                #    stride = round(npnts / n_slices)
                #    nloop = list(range(stride))
                #    for w in nloop:
                #        subset_t_col_o3 = tcolo3[w::stride]
                #        subset_t_cld = tcld[w::stride]
                #        self.add_slice(i, j, subset_t_cld, subset_t_col_o3)

    
    def calc_seasonal_means(self):
        """Calculates the mean o3 mixing ratio using Gaussian weights or counts.
        This is to be applied at the end of processing to get the final data that will be saved and plotted.
        """
        self.mean_go3vmr = np.divide(self.go3vmr, self.gwgt, where=self.gcnt != 0)
        self.mean_gerr = np.divide(self.gerr, self.gwgt, where=self.gcnt != 0)
        #self.mean_gerr = self.gerr #np.sqrt((1/self.gerr), where=self.gcnt != 0)
        self.mean_gwgt = np.divide(self.gwgt, self.gcnt, where=self.gcnt != 0)
        self.mean_cld_p_range = np.divide(self.pcld_range, self.gcnt, where=self.gcnt != 0)
        self.mean_cld_p_ceil = np.divide(self.pcld_ceil, self.gcnt, where=self.gcnt != 0)
        self.mean_go3vmr[self.gcnt == 0] = np.nan
        self.mean_gerr[self.gcnt == 0] = np.nan
        self.mean_gwgt[self.gcnt == 0] = np.nan
        self.mean_cld_p_range[self.gcnt == 0] = np.nan
        self.mean_cld_p_ceil[self.gcnt == 0] = np.nan
        self.gcnt[self.gcnt == 0] = np.nan   # Watch out for this rewriting of gcnt in the future

    def print_report(self):
        """Prints a report of useful data and reasons for data loss."""
        logging.info(f'Max no. of data points in a gridsquare: {np.nanmax(self.gcnt)}')
        # Track reasons for data loss:
        logging.info(f'(1) Too few points: {self.loss_count["too_few_points"]}')
        logging.info(f'(2) Low cloud height range: {self.loss_count["low_cloud_height_range"]}')
        logging.info(f'(3) Low cloud height std dev: {self.loss_count["low_cloud_height_std"]}')
        logging.info(f'(4) Large error: {self.loss_count["large_error"]}')
        logging.info(f'(5) Significantly less than zero: {self.loss_count["sig_diff_from_zero"]}')
        logging.info(f'(6) Outlier (o3 > 200 pptv): {self.loss_count["o3_outlier"]}')
        logging.info(f'(7) Non-uniform stratosphere: {self.loss_count["non_uni_strat"]}')
        logging.info(f'(8) Successful retrievals: {self.cloud_slice_count}')
        logging.info(f'(9) Total possible points: {(sum(self.loss_count.values()) + self.cloud_slice_count)}')
        logging.info(f'Mean % points retained for cloud slicing: {np.mean(self.postfilt)}')

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

        # Save data: UT o3 VMR (go3vmr), UT o3 error (gerr), No. of data points (gcnt)
        uto3 = ncout.createVariable('uto3', np.float32, ('lon', 'lat'))
        uto3.units = 'ppbv'
        uto3.long_name = f'O3 mixing ratios at {self.pmin}-{self.pmax} hPa obtained with cloud-slicing'
        uto3[:] = self.mean_go3vmr

        uto3err = ncout.createVariable('uto3err', np.float32, ('lon', 'lat'))
        uto3err.units = 'ppbv'
        uto3err.long_name = f'Standard error of O3 mixing ratios at {self.pmin}-{self.pmax} hPa obtained with cloud-slicing'
        uto3err[:] = self.mean_gerr

        utcldrange = ncout.createVariable('cld_top_p_range', np.float32, ('lon', 'lat'))
        utcldrange.units = 'hPa'
        utcldrange.long_name = 'Gridded mean range in cloud top pressures used to cloud-slice TROPOMI O3'
        utcldrange[:] = self.mean_cld_p_range

        utcldceil = ncout.createVariable('cld_top_p_ceil', np.float32, ('lon', 'lat'))
        utcldceil.units = 'hPa'
        utcldceil.long_name = 'Gridded mean ceiling of cloud top pressures used to cloud-slice TROPOMI O3'
        utcldceil[:] = self.mean_cld_p_ceil

        nobs = ncout.createVariable('nobs', np.float32, ('lon', 'lat'))
        nobs.units = 'unitless'
        nobs.long_name = 'Number of observations in each gridsquare used to derive O3 mixing ratios with cloud slicing'
        nobs[:] = self.gcnt

        ncout.close()
        
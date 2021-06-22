#!/usr/bin/python

''' Cloud-slicing steps applied to a cluster of data using as input the partial O3 columns in molecules/m2 and cloud top heights in hPa. 

If successful, the output is O3 mixing ratios in ppbv. Other output is the estimated error on the O3 mixing ratio and the mean cloud top pressure (hPa) for the cluster.

If the cloud-slicing step is unsuccessful, all values are NaN and a reason the cloud-slicing failed is provided.

Use as part of a Python script:
::
    from uptrop.cloud_slice_ut_o3 import cldslice, CLOUD_SLICE_ERROR_ENUM
    # Dummy functions
    partial_columns = get_partial_cloud_columns()
    cloud_top_heights = get_cloud_top_heights()

    ratio, ratio_error, exception_number, mean_cloud_error = cldslice(partial_columns, cloud_top_heights)
    if exception_number != 0:
        print(CLOUD_SLICE_ERROR_ENUM[exception_number])
    print("Cloud ratio: {}".format(ratio))
    print("Cloud ratio error: {}".format(ratio_error))
    print("Mean cloud error: {}".format(mean_cloud_error))
'''

# Import relevant packages:
import sys
import os
import numpy as np

# Import hack
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..'))

from uptrop.bootstrap import rma
from uptrop.constants import AVOGADRO as na
from uptrop.constants import G as g
from uptrop.constants import MW_AIR as mmair

from scipy import stats

CLOUD_SLICE_ERROR_ENUM = {
    1: "too_few_points",
    2: "low_cloud_height_range",
    3: "low_cloud_height_std",
    4: "large_error",
    5: "much_less_than_zero",
    6: "o3_outlier",
    7: "non_uni_strat"
}


def cldslice(pcolo3,cldtophgt):
    """
    Compute upper troposphere O3 using partial columns above
    cloudy scenes.

    Determine O3 mixing ratio by regressing O3 partial columns
    against cloud-top heights over cloudy scenes.

    :param pcolo3: vectors of partial columns in molec/m2
    :type pcolo3: list of floats
    :param cldtophgt: corresponding cloud top heights in hPa.
    :type cldtophgt: list of floats

    :return: O3 volumetric mixing ratio, corresponding estimated error on the
            cloud-sliced O3 value, a number to identify which filtering
            criteria led to loss of data in the case that the cloud-sliced
            O3 value ia nan, and the mean cloud pressure of data retained
            after 10th and 90th percentile filtering.
    :rtype: tuple
    """
    # Initialize:
    utmro3=0.0
    utmro3err=0.0
    error_state=0

    # Define factor to convert slope of O3 partial column vs pressure
    # to VMR:
    den2mr=np.divide((np.multiply(g,mmair)),na)

    # Get 10th and 90th percentiles of data population:
    #p10=np.percentile(pcolo3,10)
    #p90=np.percentile(pcolo3,90)

    # Remove outliers determined as falling outside the 10th and 90th 
    # percentile range. Not include this or instead using 5th and 95th leads
    # to overestimate in cloud-sliced UT O3 compared to the "truth":
    #sind=np.where((pcolo3>p10)&(pcolo3<p90))[0]
    # Trim the data to remove ouliers:
    #pcolo3=pcolo3[sind]
    #cldtophgt=cldtophgt[sind]

    # Cloud pressure mean:
    mean_cld_pres=np.mean(cldtophgt)

    # Get number of points in vector:
    npoints=len(cldtophgt)

    # Only consider data with more than 5 points for reasonably
    # robust statistics. This step is added to account for data loss
    # removing outliers:
    #if npoints<=10:
    #    error_state=1
    #    utmro3=np.nan
    #    utmro3err=np.nan
    #    return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Get cloud top height standard deviation:
    stdcld=np.std(cldtophgt)
    # Get cloud top height range:
    diffcld=np.nanmax(cldtophgt)-np.nanmin(cldtophgt)

    # Only consider scenes with a dynamic range of clouds:
    # (i) Cloud range:
    if diffcld<=140:
        error_state=2
        utmro3=np.nan
        utmro3err=np.nan
        return (utmro3, utmro3err, error_state, mean_cld_pres)

    # (ii) Cloud standard deviation:
    if stdcld<=30:
        error_state=3
        utmro3=np.nan
        utmro3err=np.nan
        return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Get regression statistics:
    # Partial O3 column (molec/m2) vs cloud top height (hPa):
    # 300 iterations of regression chosen to compromise between
    # statistics and computational efficiency:
    result=rma(cldtophgt*1e2,pcolo3,npoints,300)

    #if ts_slope>=0:
    #    print('RMA slope (e15): ',result[0]*1e-15)
    #    print('RMA slope rel err: ',result[2]/result[0])
    #    print('TS slope (e15): ',ts_slope*1e-15)
    #    print('TS slope rel err: ',ts_slope_err/ts_slope)
    #    print('TS slope err (e15): ',ts_slope_err*1e-15)
    #    print(ts_slope_err*1e-15)
    #    print('Slope low high: ', ts_reg[2:])
    #    r = stats.pearsonr(cldtophgt*1e2,pcolo3)
    #    print('R = ',r)
    #    if r[0]>0.6: sys.exit()

    # Calculate linear correlation:
    #corr = stats.pearsonr(cldtophgt*1e2,pcolo3)
    #    print('R = ',r)
    #if corr[0]>=0.5:
    #    # Theil-Sen regression:
    #    ts_reg = stats.theilslopes(pcolo3, cldtophgt*1e2, 0.9)
    #    ts_slope = ts_reg[0]
    #    ts_slope_err = np.abs(ts_reg[0]-ts_reg[2:])
    #else:
    #    error_state=5   # Need to update.
    #    utmro3=np.nan
    #    utmro3err=np.nan
    #    return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Remove data with relative error > 100%:
    if np.absolute(np.divide(result[2], result[0]))>1.0:
    #if np.min(np.abs(np.divide(ts_slope_err, ts_slope)))>1.0:
        error_state=4
        utmro3=np.nan
        utmro3err=np.nan
        return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Account for negative values:
    # Set points with sum of slope and error less than zero to nan.
    # This is to account for noise in the data that hover near zero.
    # This could be reduced to a single line. Eventually update to:
    #if result[0]<0 and np.add(result[0],result[2])<0):
    #    error_state=5
    #    utmro3=np.nan
    #    utmro3err=np.nan
    #    return (utmro3, utmro3err, error_state, mean_cld_pres)
    if result[0]<0 and (not np.isnan(utmro3)):
        if (np.add(result[0],result[2])<0):
    #if ts_slope<0 and (not np.isnan(utmro3)):
    #    if np.max(np.add(ts_slope,ts_slope_err))<0:
            error_state=5
            utmro3=np.nan
            utmro3err=np.nan
            return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Proceed with estimating O3 mixing ratios for retained data:
    #if not np.isnan(utmro3):
    slope=result[0]
    #slope=ts_slope
    #slope=np.multiply(slope,sf)
    slope_err=result[2]
    #slope_err=np.max(ts_slope_err)
    #slope_err=np.multiply(slope_err,sf)
    # Convert slope to mol/mol:
    utmro3=np.multiply(slope,den2mr)
    # Convert error to mol/mol:
    utmro3err=np.multiply(slope_err,den2mr)
    # Convert UT O3 from mol/mol to ppb:
    utmro3=np.multiply(utmro3,1e+9)
    # Convert UT O3 error from mol/mol to ppb
    utmro3err=np.multiply(utmro3err,1e+9)

    # Finally, remove outliers in the cloud-sliced O3
    # 200 ppbv threshold is chosen, as far from likely range.
    # Scale factor applied to TROPOMI UT O3 to account for
    # positive bias in free tropospheric O3:
    if utmro3>200:
        error_state=6
        utmro3=np.nan
        utmro3err=np.nan
        return (utmro3, utmro3err, error_state, mean_cld_pres)
    else:
        return (utmro3, utmro3err, error_state, mean_cld_pres)

#!/usr/bin/python

''' Cloud-slicing steps applied to a cluster of data using as input the partial NO2 columns in molecules/m2 and cloud top heights in hPa. 

If successful, the output is NO2 mixing ratios in pptv. Other output is the estimated error on the NO2 mixing ratio and the mean cloud top pressure (hPa) for the cluster.

If the cloud-slicing step is unsuccessful, all values are NaN and a reason the cloud-slicing failed is provided.

Use as part of a Python script:
::
    from uptrop.cloud_slice_ut_no2 import cldslice, CLOUD_SLICE_ERROR_ENUM
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
    5: "sig_diff_from_zero",
    6: "no2_outlier",
    7: "non_uni_strat"
}


def cldslice(pcolno2,cldtophgt,cld_diff_thold):
    """
    Compute upper troposphere NO2 using partial columns above
    cloudy scenes.

    Determine NO2 mixing ratio by regressing NO2 partial columns
    against cloud-top heights over cloudy scenes.

    :param pcolno2: vectors of partial columns in molec/m2
    :type pcolno2: list of floats
    :param cldtophgt: corresponding cloud top heights in hPa.
    :type cldtophgt: list of floats

    :return: NO2 volumetric mixing ratio, corresponding estimated error on the
            cloud-sliced NO2 value, a number to identify which filtering
            criteria led to loss of data in the case that the cloud-sliced
            NO2 value ia nan, and the mean cloud pressure of data retained
            after 10th and 90th percentile filtering.
    :rtype: tuple
    """
    # Initialize:
    utmrno2=0.0
    utmrno2err=0.0
    error_state=0

    # Define factor to convert slope of NO2 partial column vs pressure
    # to VMR:
    den2mr=np.divide((np.multiply(g,mmair)),na)

    # Get 10th and 90th percentiles of data population:
    #p10=np.percentile(pcolno2,10)
    #p90=np.percentile(pcolno2,90)

    # Remove outliers determined as falling outside the 10th and 90th 
    # percentile range. Not include this or instead using 5th and 95th leads
    # to overestimate in cloud-sliced UT NO2 compared to the "truth":
    #sind=np.where((pcolno2>p10)&(pcolno2<p90))[0]
    # Trim the data to remove ouliers:
    #pcolno2=pcolno2[sind]
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
    #    utmrno2=np.nan
    #    utmrno2err=np.nan
    #    return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Get cloud top height standard deviation:
    stdcld=np.std(cldtophgt)
    # Get cloud top height range:
    diffcld=np.nanmax(cldtophgt)-np.nanmin(cldtophgt)

    # Only consider scenes with a dynamic range of clouds:
    # (i) Cloud range:
    if diffcld<=cld_diff_thold:
        error_state=2
        utmrno2=np.nan
        utmrno2err=np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # (ii) Cloud standard deviation:
    if stdcld<=30:
        error_state=3
        utmrno2=np.nan
        utmrno2err=np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Get regression statistics:
    # Partial NO2 column (molec/m2) vs cloud top height (hPa):
    # 300 iterations of regression chosen to compromise between
    # statistics and computational efficiency:
    result=rma(cldtophgt*1e2,pcolno2,npoints,300)
    
    # Try Theil-Sen regressor instead:
    # Test with GEOS-Chem.    
    ts_reg = stats.mstats.theilslopes(pcolno2, cldtophgt*1e2, 0.9)
    ts_slope = ts_reg[0] 
    ts_slope_err = max(np.abs(ts_reg[0]-ts_reg[2:4])) # Conservative approach: choose largest error value

    # Compare the two:
    #print('RMA slope: ', result[0])
    #print('Theil Slope: ', ts_slope)
    #print('Theil Slope CI: ', ts_reg[2:4])
    #sys.exit()

    #if ts_slope>=0:
    #    print('RMA slope (e12): ',result[0]*1e-12)
    #    print('RMA slope rel err: ',result[2]/result[0])
    #    print('TS slope (e12): ',ts_slope*1e-12)
    #    print('TS slope rel err: ',ts_slope_err/ts_slope)
    #    print('TS slope err (e12): ',ts_slope_err*1e-12)
    #    print(ts_slope_err*1e-12)
    #    print('Slope low high: ', ts_reg[2:])
    #    r = stats.pearsonr(cldtophgt*1e2,pcolno2)
    #    print('R = ',r)
    #    if r[0]>0.6: sys.exit()

    # Remove data with relative error > factor of 3:
    #if np.absolute(np.divide(result[2], result[0]))>1.0:
    if np.divide(ts_slope_err, ts_slope)>2.0:
        error_state=4
        utmrno2=np.nan
        utmrno2err=np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Account for negative values:
    # Set points with sum of slope and error less than zero to nan.
    # This is to account for noise in the data that hover near zero.
    # This could be reduced to a single line. Eventually update to:
    #if result[0]<0 and np.add(result[0],result[2])<0):
    #    error_state=5
    #    utmrno2=np.nan
    #    utmrno2err=np.nan
    #    return (utmrno2, utmrno2err, error_state, mean_cld_pres)
    #if result[0]<0 and (not np.isnan(utmrno2)):
    #    if (np.add(result[0],result[2])<0):
    if ts_slope<0:     # and (not np.isnan(utmrno2)):
        if np.add(ts_slope,ts_slope_err)<0:
            error_state=5
            utmrno2=np.nan
            utmrno2err=np.nan
            return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Proceed with estimating NO2 mixing ratios for retained data:
    #if not np.isnan(utmrno2):
    #slope=result[0]
    slope=ts_slope
    #slope=np.multiply(slope,sf)
    #slope_err=result[2]
    slope_err=np.max(ts_slope_err)
    #slope_err=np.multiply(slope_err,sf)
    # Convert slope to mol/mol:
    utmrno2=np.multiply(slope,den2mr)
    # Convert error to mol/mol:
    utmrno2err=np.multiply(slope_err,den2mr)
    # Convert UT NO2 from mol/mol to ppt:
    utmrno2=np.multiply(utmrno2,1e+12)
    # Convert UT NO2 error from mol/mol to ppt
    utmrno2err=np.multiply(utmrno2err,1e+12)

    # Finally, remove outliers in the cloud-sliced NO2
    # 200 pptv threshold is chosen, as far from likely range.
    # Scale factor applied to TROPOMI UT NO2 to account for
    # positive bias in free tropospheric NO2:
    if np.abs(utmrno2)>200:
        error_state=6
        utmrno2=np.nan
        utmrno2err=np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)
    else:
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

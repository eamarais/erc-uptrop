''' Cloud-slicing steps applied to a cluster of data using as input the partial NO2 columns in molecules/m2 and cloud top heights in hPa. 

If successful, the output is NO2 mixing ratios in pptv. Other output is the estimated error on the NO2 mixing ratio and the mean cloud top pressure (hPa) for the cluster.

If the cloud-slicing step is unsuccessful, all values are NaN and a reason the cloud-slicing failed is provided.
'''

# Import relevant packages and modules:
import numpy as np
from scipy import stats

from .constants import (
    AVOGADRO as na,
    G as g,
    MW_AIR as mmair,
)


def cloud_slice_no2(pcolno2, cldtophgt, cld_diff_thold, min_npoints = 10, ci = 0.9,
                    relative_thold = None, check_negative = True, outlier_thold = None):
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
    utmrno2 = 0.0
    utmrno2err = 0.0
    error_state = 0

    # Define factor to convert slope of NO2 partial column vs pressure
    # to VMR:
    den2mr = np.divide((np.multiply(g, mmair)), na)

    # Cloud pressure mean:
    mean_cld_pres = np.mean(cldtophgt)

    # Get number of points in vector:
    npoints = len(cldtophgt)

    # Only consider data with more than a minimum number of points for reasonably
    # robust statistics. This step is added to account for data loss
    # removing outliers:
    if npoints <= min_npoints:
        error_state = 1
        utmrno2 = np.nan
        utmrno2err = np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Get cloud top height standard deviation:
    stdcld = np.std(cldtophgt)
    # Get cloud top height range:
    diffcld = np.nanmax(cldtophgt) - np.nanmin(cldtophgt)

    # Only consider scenes with a dynamic range of clouds:
    # (i) Cloud range:
    if diffcld <= cld_diff_thold:
        error_state = 2
        utmrno2 = np.nan
        utmrno2err = np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # (ii) Cloud standard deviation:
    if stdcld <= 30:
        error_state = 3
        utmrno2 = np.nan
        utmrno2err = np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Get regression statistics:
    # Using the input confidence interval:
    result = stats.theilslopes(pcolno2, cldtophgt * 1e2, ci)
    
    # Get the slope:
    slope = result[0]
    
    # Calculate slope error as average of errors on either side of slope:
    slope_err = 0.5 * ((result[3] - result[0]) + (result[0] - result[2]))

    # Remove data with relative error (for example, relative_thold = 100%):
    if relative_thold is not None and np.absolute(np.divide(slope_err, slope)) > relative_thold:
        error_state = 4
        utmrno2 = np.nan
        utmrno2err = np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Account for negative values:
    # Set points with sum of slope and error less than zero to nan.
    # This is to account for noise in the data that hover near zero.
    if check_negative and slope < 0:
        error_state = 5
        utmrno2 = np.nan
        utmrno2err = np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Proceed with estimating NO2 mixing ratios for retained data:
    # Convert slope to mol/mol:
    utmrno2 = np.multiply(slope, den2mr)
    # Convert UT NO2 from mol/mol to ppt:
    utmrno2 = np.multiply(utmrno2, 1e12)
    
    # Remove outliers in the cloud-sliced NO2
    # 200 pptv threshold may be chosen, as far from likely range.
    # Scale factor applied to TROPOMI UT NO2 to account for
    # positive bias in free tropospheric NO2:
    if outlier_thold is not None and utmrno2 > outlier_thold:
        error_state = 6
        utmrno2 = np.nan
        utmrno2err = np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres) 
        
    # Convert error to mol/mol:
    utmrno2err = np.multiply(slope_err, den2mr)
    # Convert UT NO2 error from mol/mol to ppt
    utmrno2err = np.multiply(utmrno2err, 1e12)
    
    # Output:
    return (utmrno2, utmrno2err, error_state, mean_cld_pres) 
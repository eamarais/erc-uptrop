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
    DU_TO_MOLECULES_PER_CM2 as du2moleccm2,
)


def cloud_slice_o3(pcolo3, cldtophgt, cld_diff_thold, pmin, pmax):
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
    # No longer apply this, as Bex Horner's comparison of cloud-sliced
    # O3 with aircraft measurements from NASA DC8 campaigns shows that
    # this causes disagreement between the two datasets:
    #p10=np.percentile(pcolo3,10)
    #p90=np.percentile(pcolo3,90)

    # Remove outliers determined as falling outside the 10th and 90th 
    # percentile range. Not include this or instead using 5th and 95th leads
    # to overestimate in cloud-sliced UT O3 compared to the "truth":
    #sind=np.where((pcolo3>p10)&(pcolo3<p90))[0]
    #sind=np.where((pcolo3<p99))[0]
    # Trim the data to remove outliers:
    #pcolo3=pcolo3[sind]
    #cldtophgt=cldtophgt[sind]

    # Cloud pressure mean:
    mean_cld_pres=np.mean(cldtophgt)

    # Get number of points in vector:
    npoints=len(cldtophgt)

    # Only consider data with more than 5 points for reasonably
    # robust statistics. This step is added to account for data loss
    # removing outliers:
    if npoints<=10:
        error_state=1
        utmro3=np.nan
        utmro3err=np.nan
        return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Get cloud top height standard deviation:
    stdcld=np.std(cldtophgt)
    # Get cloud top height range:
    diffcld=np.nanmax(cldtophgt)-np.nanmin(cldtophgt)

    # Only consider scenes with a dynamic range of clouds:
    # (i) Cloud range:
    if diffcld<=cld_diff_thold:   #140:
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
    #result=rma(cldtophgt*1e2,pcolo3,len(pcolo3),1000)

    # Using 90% confidence interval:
    result = stats.theilslopes(pcolo3, cldtophgt)#, 0.9)

    # Check for stratospheric influence following criterion in
    # Ziemke et al. (2003) that absolute total column O3 gradient is > 50 DU.
    # Ziemke warn that it's a subjective value.
    gradient = np.abs(result[0]*pmax - result[0]*pmin)
    if gradient>=50:
        error_state = 7
        utmro3=np.nan
        utmro3err=np.nan
        return (utmro3, utmro3err, error_state, mean_cld_pres)

    # calculate slope error as average of errors on either side of slope:
    #slope_err = 0.5*((result[3] - result[0])+(result[0] - result[2]))
    slope_stderr = (result[3]-result[2]) / 3.95   # 3.29 @ 90% CI; 3.92 @ 95% CI
    slope_sdev = np.sqrt(npoints) * slope_stderr

    # Pearson's correlation coefficient:
    corr = stats.pearsonr(pcolo3, cldtophgt)

    # Plot the data:
    #if (result[0]>0 and corr[0]>0.4):
        #plt.plot(cldtophgt, pcolo3, marker='o', linestyle="None")
        #plt.ylim([270,285])
        #plt.xlim([310,450])

        #print('O3: ', result[0]*1.27*1e3)
        #print('O3 error: ', slope_err*1.27*1e3)
        #print('R: ', corr[0])

        #plt.show()

    #    sys.exit()

    # Remove data with relative error > 100%:
    #if np.absolute(np.divide(slope_err, result[0]))>2.0:
    #    error_state=4
    #    utmro3=np.nan
    #    utmro3err=np.nan
    #    return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Define slope and slope error:
    slope=result[0]
    slope_err = slope_stderr

    # Account for negative values:
    # Set points with sum of slope and error less than zero to nan.
    # This is to account for noise in the data that hover near zero.
    # This could be reduced to a single line. Eventually update to:
    if (np.add(result[0],slope_err)<0):
    #if result[0]<0:
        error_state=5
        utmro3=np.nan
        utmro3err=np.nan
        return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Weak linear correlation. Use "Large error" error message for this for now:
    #if corr[0]<0.4:
    #    error_state=4
    #    utmro3=np.nan
    #    utmro3err=np.nan
    #    return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Proceed with estimating O3 mixing ratios for retained data:

    # Ziemke et al. papers 2-sigma rule:
    if np.abs(slope_err*2) >= np.abs(slope):
        error_state=4
        utmro3=np.nan
        utmro3err=np.nan
        return (utmro3, utmro3err, error_state, mean_cld_pres)

    #print('Slope: ', slope)
    #print('CI: ', result[3], result[2])
    #print('R: ', corr[0])
    #print('Points: ', npoints)
    #print('O3 vmr: ', slope*1.27*1e3)
    #print('CI: ', result[3]-result[2])
    #print('Standard error: ', slope_stderr)
    #print('Standard deviation: ', slope_sdev)
    #print('O3 vmr error: ', slope_err*1.27*1e3)
    #print(den2mr * du2moleccm2 * 1e4 * 1e-2 * 1e9)
    #print('O3 vmr (ppbv): ', slope*den2mr * du2moleccm2 * 1e4 * 1e-2 * 1e9)
    #print('O3 vmr error: ', slope_err*den2mr*du2moleccm2 * 1e4 * 1e-2 * 1e9)

    # Convert slope from DU/Pa to ppbv:
    utmro3 = slope * den2mr * du2moleccm2 * 1e4 * 1e-2 * 1e9
    #utmro3 = slope * 1.27 * 1e3
    # Convert slope error from DU/Pa to ppbv:
    utmro3err = slope_err * den2mr * du2moleccm2 * 1e4 * 1e-2 * 1e9
    #utmro3err = slope_sdev * 1.27 * 1e3
    
    # Convert slope to mol/mol:
    #utmro3=np.multiply(slope,den2mr)
    # Convert error to mol/mol:
    #utmro3err=np.multiply(slope_err,den2mr)
    # Convert UT O3 from mol/mol to ppb:
    #utmro3=np.multiply(utmro3,1e+9)
    # Convert UT O3 error from mol/mol to ppb:
    #utmro3err=np.multiply(utmro3err,1e+9)

    # Put in 2-sigma condition:
    # if 2*utmro3err > utmr03, then set values as np.nan.

    # Output:
    return (utmro3, utmro3err, error_state, mean_cld_pres)

    # Finally, remove outliers in the cloud-sliced O3
    # 200 pptv threshold is chosen, as far from likely range.
    # Scale factor applied to TROPOMI UT O3 to account for
    # positive bias in free tropospheric O3:
    #if utmro3>200:
    #    error_state=6
    #    utmro3=np.nan
    #    utmro3err=np.nan
    #    return (utmro3, utmro3err, error_state, mean_cld_pres)
    #else:
Overview of Capabilities
========================

Uptrop-Py is python coded software used to apply the so-called cloud-slicing technique to satellite observations of column densities of nitrogen dioxide (NO\ :sub:`2`) to retrieve vertical profiles of mixing ratios of NO\ :sub:`2` throughout the troposphere on a global scale. 

Code output includes data in NetCDF format, plots in postscript format, and metrics at the end of the log file. Sample plots and logfile metrics are in the `Gallery <https://erc-uptrop.readthedocs.io/en/latest/Gallery.html>`__.

So far the code has been applied to Sentinel-5P TROPOMI instrument obesrvations of total columns of NO\ :sub:`2` using cloud information from two distinct TROPOMI cloud products to obtain NO\ :sub:`2` concentrations in a single layer in the upper troposphere from 450 hPa (~8 km) to 180 hPa (~12 km). 

The steps in the retrieval algorithm are detailed in the open access peer-reviewed `Marais et al. (2021) <https://doi.org/10.5194/amt-14-2389-2021>`__ paper published in Atmospheric Measurement Techniques. 

The data generated for this paper are available for public download from the `UCL Data Repository <https://doi.org/10.5522/04/14586558.v1>`__.

A static version of the software that was used in the above paper to generate the above data is hosted on Zenodo and available for public download and use as `version 1.1.0 <https://doi.org/10.5281/zenodo.4058442>`__. 

About Uptrop
============

**Uptrop** is a Python-based software designed to apply the cloud-slicing technique to
satellite observations of column densities of nitrogen dioxide (NO\ :sub:`2`) and
ozone (O\ :sub:`3`) to retrieve vertical profiles of mixing ratios of NO\ :sub:`2`
and O\ :sub:`3` throughout the troposphere on a global scale.

The software was developed using total column densities of NO\ :sub:`2` from the
ESA Sentinel-5P TROPOMI instrument to derive NO\ :sub:`2` mixing ratios in the global upper
troposphere (8-12 km) at ~100 km spatial resolution. The steps in the retrieval algorithm
are detailed in the open access peer-reviewed 
`Marais et al. (2021) <https://doi.org/10.5194/amt-14-2389-2021>`__ paper published
in Atmospheric Measurement Techniques. The data generated for this paper are available
for public download from the `UCL Data Repository <https://doi.org/10.5522/04/14586558.v1>`__.
A static version of the software that was used in the above paper to generate the
above data is hosted on Zenodo and available for public download and use as
`version 1.1.0 <https://doi.org/10.5281/zenodo.4058442>`__.

The software is now updated to retrieve NO\ :sub:`2` and O\ :sub:`3` mixing ratios in
custom vertical layers throughout the troposphere using Sentinel-5P TROPOMI instrument
obesrvations of total columns of NO\ :sub:`2` and O\ :sub:`3`, and cloud information from two
distinct TROPOMI cloud products.

Code output includes data in NetCDF format, plots in postscript format, and metrics
at the end of the log file. Sample plots and logfile metrics are in the
`Gallery <https://erc-uptrop.readthedocs.io/en/latest/Gallery.html>`__.
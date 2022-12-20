What is Uptrop-Py?
==================

A library for applying the cloud-slicing technique to TROPOMI total
columns to obtain upper tropospheric (450-180 hPa) mixing ratios of NO2.
This includes a feasibility test of the technique using synthetic
observations from GEOS-Chem, validation and quantification of bias
corrections by comparing TROPOMI and Pandora total columns at
high-altitide sites, application of cloud-slicing using either the
FRESCO-S or ROCINN-CAL cloud products, and comparison of effective cloud
fractions and cloud top pressures from the two cloud products.

In the modules, the FRESCO-S product is referred to as FRESCO and the
ROCINN-CAL product as DLR.

A description of the data used, the sources of the data, the details of
the GEOS-Chem model simulations, and other relevant information for
developing the codes can be found in Marais et al. (2020).

The easiest way to use the library is to call the scripts from a
terminal; see the rest of this page for examples. You can also import
the components into your own scripts - see the individual module
descriptions linked in the sidebar.

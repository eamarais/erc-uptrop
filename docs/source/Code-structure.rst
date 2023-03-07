Code Structure
================

===========
Core Code:
===========

Used for cloud-slicing satellite observations of total column densities of NO2:

.. option:: ./erc-uptrop/uptrop/cloud_slice_tropomi_no2.py

Module that reads in TROPOMI data, calls the cloud slicing routine, grids data, calculates seasonal means, outputs data and sample plots.

.. option:: ./erc-uptrop/uptrop/cloud_slice_no2.py

Routine that cloud-slices the TROPOMI NO2 data prepared in :file:`cloud_slice_tropomi_no2.py`. Calculates NO2 concentration and associated error in pptv.

.. option:: ./erc-uptrop/uptrop/height_pressure_converter.py

Routine to convert cloud top height to cloud top pressure. Used for the TROPOMI CLOUD_OFFL product only.

.. option:: ./erc-uptrop/uptrop/boostrap.py

Routine that calculates the reduced major axis regression slope and intercept values and errors estimated with bootstrapping. Adapted into Python from the IDL `GAMAP package <http://wiki.seas.harvard.edu/geos-chem/index.php/General_GAMAP_usage>`__.

.. option:: ./erc-uptrop/uptrop/gamap_colormap.py

The white-yellow-green-orange-red colorbar from the `GCPy Python package <https://gcpy.readthedocs.io/en/stable/>`__.

.. option:: ./erc-uptrop/uptrop/constants.py

Constants and conversion factors used for cloud-slicing. 

===============
Additional Code:
===============

Used by `Marais et al. (2021) <https://amt.copernicus.org/articles/14/2389/2021/>`__ to assess the TROPOMI NO2 columns and cloud products:

.. option:: ./erc-uptrop/uptrop/fresco_cld_error.py

Routine to compare two TROPOMI cloud products used to cloud-slice TROPOMI NO2. 

.. option:: ./erc-uptrop/uptrop/read_pandora.py

Routine to read in ground-based Pandora total and tropospheric NO2 column data in the format provided by the `Pandonia Global Network <https://www.pandonia-global-network.org/>`__. 

.. option:: ./erc-uptrop/uptrop/compare_tropomi_pandora.py

Routine to sample and compare coincident data from TROPOMI and Pandora total and tropospheric NO2 columns at high-altitude sites.

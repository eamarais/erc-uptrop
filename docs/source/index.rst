.. Uptrop documentation master file, created by
   sphinx-quickstart on Wed Jul  1 13:29:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Uptrop's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   bootstrap
   cloud_slice_ut_no2
   compare_tropomi_pandora
   constants
   convert_height_to_press
   fresco_cld_err
   gamap_colormap
   read_pandora
   tropomi_ut_no2
   ut_no2_gc_test



Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
============
A library for applying the cloud-slicing technique to TROPOMI total columns to obtain upper tropospheric (450-180 hPa) mixing ratios of NO2. This includes a feasibility test of the technique using synthetic observations from GEOS-Chem, validation and quantification of bias corrections by comparing TROPOMI and Pandora total columns at high-altitide sites, application of cloud-slicing using either the FRESCO-S or ROCINN-CAL cloud products, and comparison of effective cloud fractions and cloud top pressures from the two cloud products.

Requirements
------------
- Python >= 3.6
- Numpy
- NetCDF4
- Dateutil
- Matplotlib
- Basemap


If you are using Conda, an environment.yml file is provided.

Example of use
--------------

.. code-block:: bash

   python tropomi_ut_no2 /path/to/tropomi/data/folder/ /path/to/output/folder

For more calling details


Scripts
-------

``ut_no2_gc_test``

.. automodule:: uptrop.ut_no2_gc_test

``cloud_slice_ut_no2``

.. automodule:: uptrop.cloud_slice_ut_no2

``compare_tropomi_pandora``

.. automodule:: uptrop.compare_tropomi_pandora

``fresco_cld_err``

.. automodule:: uptrop.fresco_cld_err

``tropomi_ut_no2``

.. automodule:: uptrop.tropomi_ut_no2
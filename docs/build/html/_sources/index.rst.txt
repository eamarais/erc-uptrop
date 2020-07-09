.. Uptrop documentation master file, created by
   sphinx-quickstart on Wed Jul  1 13:29:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Uptrop's documentation!
==================================
Here you will find all the relevant details of the functionality of the code used to generate results described and discussed Marais et al. (2020).

The python code was first developed by Eloise Marais and substantially improved by John Roberts. Roberts also constucted the code documentation. 

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

In the modules, the FRESCO-S product is referred to as FRESCO and the ROCINN-CAL product as DLR.

A description of the data used, the sources of the data, the details of the GEOS-Chem model simulations, and other relevant information for developing the codes can be found in Marais et al. (2020).

Requirements
------------
- Python >= 3.6
- Numpy
- NetCDF4
- Dateutil
- Matplotlib
- Basemap


If you are using Conda, an environment.yml file is provided.

Examples of use
---------------

.. code-block:: bash

   To cloud-slice synthetic partial columns from GEOS-Chem native North America (NA) nested domain and obtain seasonal mean NO2 mixing ratios at 4x5:
   python ut_no2_gc_test.py --out_path="/path/to/output/file" --resolution="4x5" --region="NA"
   
   To co-sample daily Pandora and TROPOMI total columns of NO2 at the Izana high-altitude site:
   python compare_tropomi_pandora.py /path/to/tropomi/data/folder/ /path/to/pandora/data/folder/ /path/to/output/folder/ --no2_col="Tot" --pandora_site="izana"
   
   To cloud-slice TROPOMI partial columns of NO2 using the FRESCO-S cloud product "fresco" for June-August ("jja"):
   python tropomi_ut_no2 /path/to/tropomi/data/folder/ /path/to/output/folder/ --cloud_product="fresco" --season="jja"
   
   To obtain coincident cloud top pressures and cloud fractions from FRESCO-S and ROCINN-CAL as monthly averages for January at 1x1 using cloud top height rather than cloud top pressure from the ROCINN-CAL product:
   python fresco_cld_err.py --month="01" --number_of_days=31 --out_res="1x1" --dlr_cld_top="height"  


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

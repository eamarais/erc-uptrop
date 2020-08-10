Welcome to Uptrop's documentation!
==================================
Here you will find all the relevant details of the functionality of the code used to generate results described and discussed Marais et al. (2020).

The python code was developed for the ERC Starting Grant project UpTrop (https://maraisresearchgroup.co.uk/uptrop.html) led by Eloise Marais. Substantive code edits and improvements and documentation preparation are by software engineer John Roberts. 

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

The easiest way to use the library is to call the scripts from a terminal; see the rest of this page for examples. You can also import the components into your own scripts - see the individual module descriptions linked in the sidebar.

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

Cloud-slice synthetic partial columns from GEOS-Chem over North America at 4x5 in Spring:

.. code-block:: bash

   python ut_no2_gc_test.py --gc_dir "path/to/geoschem/data"  --out_path "/path/to/output/file" --resolution 4x5 --season mam --region NA
   
Co-sample daily Pandora and TROPOMI total NO2 columns at Izana in the first week of May:

.. code-block:: bash

   python compare_tropomi_pandora.py --trop_dir "/path/to/tropomi/data" --pan_dir "/path/to/pandora/data" --out_dir "/path/to/output/dir" --start_date 01-05-2019 --end_date 07-05-2019 --apply_bias_correction True --no2_col Trop --pandora_site izana
   
Cloud-slice Jun-Aug TROPOMI partial NO2 columns using the FRESCO-S cloud product in Summer:

.. code-block:: bash

   python tropomi_ut_no2.py --trop_dir "/path/to/tropomi/data/folder/" --out_dir "/path/to/output/folder/" --cloud_product fresco --season jja
   
Obtain January mean FRESCO-S and ROCINN-CAL cloud top pressures and fractions at 1x1:

.. code-block:: bash

   python fresco_cld_err.py  --trop_dir "/path/to/tropomi/data/folder/" --out_dir "/path/to/output/folder/" --start_date 01-01-2019 --end_date 31-01-2019 --out_res 1x1 --dlr_cld_top height

For a more advanced example bash script that runs each script over a given date range, see ``tests/integration_test.sh`` in your Uptrop folder.

Conventions
-----------

* All scripts take a path to a directory containing the relevant data (``--trop_dir``, ``--pan_dir`` and ``--gc_dir``)

* All scripts need either a ``--season`` or a ``--start_date`` and ``--end_date`` provided. If all are given, ``--season`` will override start_date/end_date

* If ``--season`` is given, then this script runs for that season 2019-2020:

        ``djf``; dececmber-januarary-februray 2019-2020

        ``mam``; march-april-may 2020

        ``jja``; june-july-august 2020

        ``son``; september-october-november 2020

* ``--start_date`` and ``--end_date`` are given as dd-mm-yyyy; for example, 20th Januarary 2019 is 20-01-2019

Scripts
-------

ut_no2_gc_test
^^^^^^^^^^^^^^

.. automodule:: uptrop.ut_no2_gc_test

cloud_slice_ut_no2
^^^^^^^^^^^^^^^^^^

.. automodule:: uptrop.cloud_slice_ut_no2

compare_tropomi_pandora
^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: uptrop.compare_tropomi_pandora

fresco_cld_err
^^^^^^^^^^^^^^

.. automodule:: uptrop.fresco_cld_err

tropomi_ut_no2
^^^^^^^^^^^^^^

.. automodule:: uptrop.tropomi_ut_no2

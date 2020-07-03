.. Uptrop documentation master file, created by
   sphinx-quickstart on Wed Jul  1 13:29:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Uptrop's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
------------
A library for comparing Sentinel 5P and Pandora atmospheric N02
I think. [ELOISE, CAN YOU PUT SOME SCIENCE HERE PLEASE I DONT REALLY KNOW WHAT THIS PROGRAM DOES]

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
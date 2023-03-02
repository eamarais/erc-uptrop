Usage
==================

Instructions below assume python is being simulated from the command line using a conda virtual environment and python version 3.9.10.

Dependencies are listed in the environment.yml file that is downloaded with the source code. 
These include:
   Numpy
   NetCDF4
   Dateutil
   Matplotlib
   Basemap

To cloud-slice TROPOMI NO2 to obtain seasonal mean NO2 mixing ratios for June-August 2019 at 180-320 hPa on a 1 degree x 1 degree global grid,
enter the following at the command line:

.. code-block:: console
   $ python cloud_slice_tropomi_no2.py --trop_dir="path-to-tropomi-data/" --out_dir="path-to-output-directory/" --cloud_product=fresco-wide --no2_prod="OFFL" --cloud_threshold="07" --grid_res="1x1" --year="2019" --pmax="180" --pmin="450" --season="jja" > log_file
   
All arguments are input as strings. 

Default input arguments are:


Options include:

FResolution options are 1x1, 2x2.5, 4x5 (degrees latitude by degrees longitude)

TROPOMI satellite data being read in assumed to be in the format ./NO2_$no2_prod/$year/$month/, where year and month are extracted from the season input information.
Output directory assumes to have ./Data/ subdirectory to place NetCDF files and ./Images/ subdirectory to place plots. 


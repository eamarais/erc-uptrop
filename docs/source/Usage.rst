Usage
==================

Instructions below assume python is being simulated from the command line using a conda virtual environment and python version 3.9.10.

Dependencies are listed in the environment.yml file that is downloaded with the source code. 
These include:
   Numpy,
   NetCDF4,
   Dateutil,
   Matplotlib,
   Basemap

To cloud-slice TROPOMI NO2 to obtain seasonal mean NO2 mixing ratios for June-August 2019 at 180-320 hPa on a 1 degree x 1 degree global grid,
enter the following at the command line:

.. code-block:: console
   $ python cloud_slice_tropomi_no2.py --trop_dir="path-to-tropomi-data/" --out_dir="path-to-output-directory/" --cloud_product="fresco-wide" --no2_prod="OFFL" --cloud_threshold="07" --grid_res="1x1" --year="2019" --pmax="180" --pmin="450" --season="jja" > log_file
   
All arguments are input as strings. 

Default input arguments are:
  --grid_res="1x1", --cloud-product="fresco-wide", --cloud-threshold="07", --pmin="180", --pmax="450", --no2-prod="PAL", --version="v1"

Resolution options are limited to:
   1x1, 2x2.5, 4x5 (degrees latitude by degrees longitude).

The TROPOMI satellite data being read in assumed to be stored in the format ./NO2_$no2_prod/$year/$month/, where year and month are extracted from the season input information.

Output directory assumes to have ./Data/ subdirectory to store NetCDF files and ./Images/ subdirectory to store plots. 

The --cloud_threshold input argument should be "07", "08", "09", "10" to limit use to optically thick clouds and associated contamination from the target compound below clouds.

Recommended cloud top pressure ranges to use include 180-320 hPa, 320-450 hPa, 450-600 hPa, 600-800 hPa, and 800-1000 hPa. The code also offers the flexibility for the user to decide on cloud top height ranges difference to those recommended, but the difference between --pmin and --pmax must exceed 100 hPa. If not, the code will stop with an error message. 

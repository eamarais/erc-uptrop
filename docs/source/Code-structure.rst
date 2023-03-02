Code Structure
================

Core code for cloud-slicing:

.. option:: erc-uptrop/uptrop/cloud_slice_tropomi_no2.py

  Module to read in TROPOMI data, call the cloud slicing routine, grid the data to seasonal means on the target grid, output data and sample plots.


erc-uptrop/uptrop/cloud_slice_tropomi_no2.py
cloud_slice_no2.py
height_pressure_converter.py
boostrap.py
gamap_colormap.py
constants.py

Additional code used in Marais et al. (2021):

fresco_cld_error.py
read_pandora.py
compare_tropomi_pandora.py

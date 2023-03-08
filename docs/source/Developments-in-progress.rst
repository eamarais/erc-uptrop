Code Updates in Progress
===========================

A new and improved version of the software is under development. It includes many updates that are categorised below.

##############
Data / Output
##############

* Change regression fit from reduced major axis to Theil, so that cloud-slicing regression fit is less impacted by non-uniformally distributed data leading to an overestimate in mixing ratios over remote locations.

* No longer bias correct tropospheric columns of NO\ :sub:`2`, as this bias correction was due to application of the reduced major axis regression fit to non-uniformally distributed data.

##############
Increased Flexibility
##############

* Adjust the code to be able to apply it to multiple pollutants. NOt just NO2. Currenlty under development is application to TROPOMI ozone to retrieve veritcal profiles of tropospheric ozone, a potent greenhouse gas.

* Ability to obtain mixing ratios in all layers in the troposphere, with the flexibility for the user to specify the pressure range, limited to a range that is at least 100 hPa for the cloud-slicing routine to succeed.

* Enable flexibility to select any start and end date, not just limited to processing data for a single season.


##############
Structural Changes
##############

* Replace Basemap with Cartopy for generating sample plots.

* Restructure code directories

* Improve variable names in the Python code and in the output NetCDF data file


If you're interested in using a version of the software currently under development, please reach out to the developers via the `GitHub Issues <https://github.com/eamarais/erc-uptrop/issues>`__ page.

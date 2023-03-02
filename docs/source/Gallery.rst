Gallery
==================

Plots output with uptrop using cartopy include three panels of global maps showing cloud-sliced nitrogen dioxide (NO2) mixing ratios at the pressure range of interest (top), the estimated error on the cloud-sliced NO2 mixing ratios, and the number of cloud-sliced values. These are not publication quality plots, but are merely for sanity checking and benchmarking.

The example below is from cloud-slicing TROPOMI NO2 in June-August 2019 at 450-180 hPa (~8-12 km):

.. figure:: images/cloudslice-tropomi-no2-sample-plot.png
  :class: with-border
  :width: 340
  :alt: Sample of plots generated with the python script



Also output to a log file with each cloud-slicing routine are diagnostics that track the maximum number of satellite pixels in the target grid, the number of satellite pixels removed in each data filtering step, the total number of successful cloud-slicing retrievals compared to the total number that could have been retrieved, and the percent of total TROPOMI pixels used for cloud slicing:

.. code-block:: text

  Max no. of data points in a gridsquare:  64.0
  (1) Too few points:  280605
  (2) Low cloud height range:  260019
  (3) Low cloud height std dev:  2105
  (4) Large error:  0
  (5) Significantly less than zero:  15664
  (6) Outlier (NO2 > 200 pptv):  0
  (7) Non-uniform stratosphere:  133461
  (8) Successful retrievals:  83516
  (9) Total possible points:  775370
  Mean % points retained:  2.141713715255334

.. figure:: images/uptrop-logfile-output.png
  :class: with-border
  :width: 380
  :alt: Sample of end of cloud-slicing log file

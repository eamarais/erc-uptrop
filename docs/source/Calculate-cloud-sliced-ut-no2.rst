Calculate cloud sliced ut no2
==================
Cloud-slicing steps applied to a cluster of data using as input the partial NO2 columns in molecules/m2 and cloud top heights in hPa.

If successful, the output is NO2 mixing ratios in pptv. Other output is the estimated error on the NO2 mixing ratio and the mean cloud top pressure (hPa) for the cluster.

If the cloud-slicing step is unsuccessful, all values are NaN and a reason the cloud-slicing failed is provided.

Use as part of a Python script:

::

   from uptrop.cloud_slice_ut_no2 import cldslice, CLOUD_SLICE_ERROR_ENUM
   # Dummy functions
   partial_columns = get_partial_cloud_columns()
   cloud_top_heights = get_cloud_top_heights()

   ratio, ratio_error, exception_number, mean_cloud_error = cldslice(partial_columns, cloud_top_heights)
   if exception_number != 0:
       print(CLOUD_SLICE_ERROR_ENUM[exception_number])
   print("Cloud ratio: {}".format(ratio))
   print("Cloud ratio error: {}".format(ratio_error))
   print("Mean cloud error: {}".format(mean_cloud_error))

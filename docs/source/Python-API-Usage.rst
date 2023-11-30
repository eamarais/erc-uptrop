Python API Usage
================


.. code-block:: python

    """
    Python API Usage Example for Uptrop

    This script demonstrates the basic usage of the Uptrop Python API to cloud-slice satellite observations.
    """

    import uptrop as up

    up.cloud_slice_tropomi(
        species='NO2',
        grid_res='1x1',
        cloud_product='fresco-wide',
        cloud_threshold='07',
        version='v1',
    )




The example script for cloud slicing TROPOMI NO2 can be found in our GitHub repository. You can download it directly using the following link:

`Download the cloud_slice_tropomi_no2_example.py script <https://github.com/eamarais/erc-uptrop/blob/refactor/docs/examples/cloud_slice_tropomi_no2_example.py>`_

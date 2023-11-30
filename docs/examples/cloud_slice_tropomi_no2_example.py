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
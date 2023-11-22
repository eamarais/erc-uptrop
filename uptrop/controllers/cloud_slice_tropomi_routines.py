""" 
Decide the cloud slice routine for the chosen species
"""

""" 
Convert the script to a class CloudSliceTropom before release
"""

import os
import datetime as dt
import logging
from dateutil import rrule as rr
import numpy as np
from .. import models
from ..utils.dates_files_utils import (
    get_date,
    season_to_date,
    get_file_list
    )


def get_class(species, class_type):
    """
    Import the classes for the chosen species.
    """
    class_name = f"{class_type}{species.upper()}"
    return getattr(models, class_name)


def process_dates(start_date, end_date, season, year):
    if start_date and end_date:
        if season:
            raise ValueError("Please provide either start_date and end_date or season, not both.")
        start_date = dt.datetime.strptime(start_date, "%Y-%m-%d")
        end_date = dt.datetime.strptime(end_date, "%Y-%m-%d")
    elif season and not (start_date or end_date):
        start_date, end_date = season_to_date(season, np.int64(year))
    else:
        raise ValueError("Please provide either start_date and end_date or season, not both.")
    return start_date, end_date


def get_and_validate_files(trop_dir, species, product, start_date, end_date, cloud_product):
    date_range = rr.rrule(rr.DAILY, dtstart=start_date, until=end_date)
    trop_files = get_file_list(trop_dir, species, product, date_range)
    logging.info(f'Found total of {len(trop_files)} files')

    # Decide which cloud files to use
    if cloud_product == "fresco-wide":
        cloud_files = trop_files
    elif cloud_product in {"dlr-ocra", "o22cld"}:
        cloud_files = get_file_list(trop_dir, "cloud", cloud_product, date_range)
    else:
        raise ValueError("Invalid cloud product; can be fresco-wide or o22cld")

    # Check for instances where there are cloud files for a TROPOMI swath, but no NO2 files and vice versa
    if cloud_product == "dlr-ocra" and len(cloud_files) != len(trop_files):
        for i, (trop_file, cloud_file) in enumerate(zip(trop_files, cloud_files)):
            trop_date = get_date(trop_file)
            cloud_date = get_date(cloud_file)
            if cloud_date < trop_date:
                del cloud_files[i]
                i = 0  # Restart iteration

    return trop_files, cloud_files


def process_grid_res(grid_res):
    grid_res_map = {
        '05x05': (0.5, 0.5),
        '1x1': (1, 1),
        '2x25': (2, 2.5),
        '4x5': (4, 5),
    }
    if grid_res not in grid_res_map:
        raise ValueError("Invalid grid; values can be 05x05, 1x1, 2x25, 4x5")
    return grid_res_map[grid_res]


# This function needs to be refactorized
# Encapsulate some of the logic into functions
# Break down this function so that other routines can be produced easily
# DRY principle!
def cloud_slice_tropomi(
    species=None,
    trop_dir=None,
    out_dir=None,
    season=None,
    year=None,
    start_date=None,
    end_date=None,
    grid_res='1x1',
    cloud_product='fresco-wide',
    cloud_threshold='07',
    fresco_440=False,
    pmin=None,
    pmax=None,
    product=None,
    version='v1',
):
    """
    This function is planned to be provided at the package level.
    The corresponding CLI tool "cloud_slice_tropomi" which takes arguments using argparse
    is built in main.py and provided in "setup.py"
    """
    # Find tropomi files for the selected dates
    start_date, end_date = process_dates(start_date, end_date, season, year)
    trop_files, cloud_files = get_and_validate_files(trop_dir, species, product, start_date, end_date, cloud_product)
    logging.info(f'Found total of {len(trop_files)} TROPOMI files and {len(cloud_files)} Cloud files')

    # Parsing cloud threshold
    cloud_threshold = float(cloud_threshold)/10

    # Parsing bollean to use or not use FRESCO-S cloud fraction at 440 nm:
    fresco_440 = fresco_440

    # Select the classes for conduncting the retrievals for the chosen species
    TropomiData = get_class(species, "Tropomi")
    CloudData = get_class(species, "CloudDataTropomi")
    GridAggregator = get_class(species, "GridAggregatorTropomi")

    # Start the cloud slicing routine
    dellat, dellon = process_grid_res(grid_res)
    grid_aggregator = GridAggregator(dellat, dellon, pmin, pmax, first = 0)

    # Get the lowercase of species name
    species_lower = species.lower()

    # Process each TROPOMI file and its corresponding cloud file
    for trop_file, cloud_file in zip(trop_files, cloud_files):

        # Process tropomi data
        trop_data = TropomiData(trop_file, product)

        # Process cloud data
        # Currently CloudDataTropomiNO2 and CloudDataTropomiO3 have different parameters for instantiation
        # For now, process these two differently
        if species_lower == "no2":
            cloud_data = CloudData(cloud_file, fresco_440, product, data_type=cloud_product)
        elif species_lower == "o3":
            cloud_data = CloudData(cloud_file, fresco_440)

        # Check the clould product and the parity of TROPOMI and cloud data
        if cloud_product == "dlr-ocra" and not cloud_data.data_parity:
            continue

        # Retrieve vertical profiles

        # Decide calc_geo_column or not based on species
        if species_lower == "no2":
            trop_data.calc_geo_column()
        elif species_lower == "o3":
            continue

        # Conduct the rest of the retrievals
        trop_data.cloud_filter_and_preprocess(cloud_data, cloud_threshold, pmax, pmin)
        grid_aggregator.initalise_grid()
        grid_aggregator.grid_trop_data(trop_data)
        grid_aggregator.apply_cloud_slice()

    # After processing all files, calculate the mean
    grid_aggregator.calc_seasonal_means()
    
    # Notify the end of the retrieval
    grid_aggregator.print_report()

    # Define output file names:
    # Define cloud product string for file name:
    str_cld_prod = cloud_product
    if fresco_440:
        str_cld_prod = 'fresco-440nm'
        
    # Data file and plot file names  
    out_file_id = (
        f'tropomi-{species}-{str_cld_prod}-{cloud_threshold}-{grid_res}-'
        f'{start_date.strftime("%Y%m%d")}-{end_date.strftime("%Y%m%d")}-'
        f'{pmin}-{pmax}hPa-{version}'
    )
    
    # Save out data (separate plotting from these)
    out_data_file = f'{out_file_id}.nc'
    out_data_file_path = os.path.join(out_dir, 'Data', out_data_file)
    logging.info(f'Saving data to: {out_data_file}')
    grid_aggregator.save_to_netcdf(out_data_file_path)

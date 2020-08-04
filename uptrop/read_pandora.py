#!/usr/bin/python

"""A small set of functions for converting a Pandora file to a dataframe.

The main function in this module is read_pandora: ::

    from uptrop.read_pandora import read_pandora
    location, pandora_df = read_pandora("pandora_file", "Tot")
    pandora_df.plot('day', 'gc_data')

The rest are ancillary functions for reading pandora data files.
"""

# Import relevant packages:
import glob
import sys
import os
import re

import numpy as np
import pandas as pd


def get_column_description_index(file_path):
    """Returns a dictionary of {description:column index} for a pandora file

    See https://regex101.com/r/gAjFtL/1 for an in-depth explanation of the regex used

    Returns two groups; group 1 is the column number, group 2 is the description.

    :param file_path: The path to the pandora data
    :type file_path: str

    :returns: The dictionary of {description:column index} suitable for passing to get_column_from_description
    :rtype: dict{str:int}
    """
    # See https://regex101.com/r/gAjFtL/1 for an in-depth explanation of this regex
    # Returns two groups; group 1 is the column number, group 2 is the description.
    searcher = re.compile(r"Column ([0-9]+): (.*)")
    with open(file_path, 'r', encoding="Latin-1") as pandora_file:
        pandora_text = pandora_file.read()
    groups = searcher.findall(pandora_text)
    column_dict = {column_description: int(column_index) for column_index, column_description in groups}
    return column_dict


def get_start_of_data(file_path):
    """Gets the line number of the start of the data itself

    Inspects the file line-by-line, counting lines until it finds the second dotted line

    :param file_path: Path to the Pandora file
    :type file_path: str

    :returns: The 1-indexed line number of the start of the data
    :rtype: int
    """
    with open(file_path, encoding="Latin-1") as pandora_file:
        # Line numbers are 1-indexed....
        line_number = 1
        # Look for the dotted line twice
        while not pandora_file.readline().startswith("-------"):
            line_number += 1
        while not pandora_file.readline().startswith("-------"):
            line_number += 1
    # Increment one more time to get the index of the first line of actual data
    line_number += 1
    return line_number


def get_column_from_description(column_dict, description):
    """Searched through the output from get_column_description_index for a given description

    :param column_dict: The output from get_column_description_index
    :type column_dict: dict
    :param description: A substring of the column description you want to find the index for

    :returns: The first index corresponding to description
    :rtype: int
    """
    index = [value for key, value in column_dict.items() if description in key][0]
    if index is None:
        return
    else:
        return index


def get_lat_lon(pandora_filepath):
    """Returns a dictionary of lat, lon extracted from the pandora file

    :param pandora_filepath: The pandora file
    :type pandora_filepath: str

    :returns: A dict of {"lat":lat, "lon":lon}
    :rtype: dict{str:int}
    """
    lat = None
    lon = None
    with open(pandora_filepath, 'r', encoding="Latin-1") as pandora_file:
        while (lat == None or lon == None):
            current_line = pandora_file.readline()
            if current_line.startswith("Location latitude [deg]:"):
                lat = float(current_line.split(":")[1])
            elif current_line.startswith("Location longitude [deg]:"):
                lon = float(current_line.split(":")[1])
            elif current_line.startswith("--------"):
                # TODO: Maybe change for exception if this might happen
                print("Lat/lon not found in file {}".format(pandora_filepath))
                return
    return {"lat": lat, "lon": lon}


def read_pandora(pandora_filepath, no2col):

    """Reads position and data from a Pandora file

    Returns two values: a dictionary of position and a dataframe of Pandora data

    Pandora data can be either all values or troposphere only

    The dictionary has key {'lat,'lon'}

    The dataframe has column headings:

    jday, sza, gc_data, no2err, qaflag, fitflag, year, month, day, hour_utc, minute


    :param pandora_filepath: The path to the pandora file
    :type pandora_filepath: str
    :param no2col: Whether to get all values 'Tot' or tropospheric values only 'Trop'
    :type no2col: str

    :returns: A tuple of the position dictionary and the dataframe
    :rtype: tuple(dict, pandas.dataframe)
    """

    loc = get_lat_lon(pandora_filepath)

    column_dict = get_column_description_index(pandora_filepath)
    dateind = get_column_from_description(column_dict, 'UT date and time for center of m')
    jdayind = get_column_from_description(column_dict, 'Fractional days since 1-Jan-2000')
    # SZA:
    szaind = get_column_from_description(column_dict,  'Solar zenith angle for center of')
    # NO2 column:
    # (a) Total:
    if no2col == 'Tot':
        no2ind = get_column_from_description(column_dict,  'Nitrogen dioxide total vertical ')
    # (b) Troposphere:
    if no2col == 'Trop':
        no2ind = get_column_from_description(column_dict,  'Nitrogen dioxide tropospheric v')
    # NO2 column error:
    # (a) Total:
    if no2col == 'Tot':
        errind = get_column_from_description(column_dict,  'Uncertainty of nitrogen dioxide total ver')
    # (b) Troposphere:
    if no2col == 'Trop':
        errind = get_column_from_description(column_dict,  'Uncertainty of nitrogen dioxide troposph')
    # Data quality flag:
    qaflagind = get_column_from_description(column_dict,  'L2 data quality flag for nitrog')

    # Level 2 fit flag:
    # There are two entries (min and max) of this in the tropospheric
    # column file. The minimum is being used here.
    fitflagind = get_column_from_description(column_dict,'Level 2 Fit data quality flag')

    data_start = get_start_of_data(pandora_filepath)

    names = ["ut_date", "jday", "sza", "no2", "no2err", "qaflag", "fitflag"]
    columns = [dateind, jdayind, szaind, no2ind, errind, qaflagind, fitflagind]
    columns = [column -1 for column in columns]  # Pandora columns are 1-indexed, Pandas are 0

    # TODO: Maybe set ut_date up as an index for easier slicing of data for the future
    df = pd.read_csv(pandora_filepath,
                     sep=" ",
                     skiprows=data_start,
                     usecols=columns,
                     names=names,
                     parse_dates=['ut_date']
                     )

    date_df = pd.DataFrame({
        "year": df.ut_date.dt.year,
        "month": df.ut_date.dt.month,
        "day": df.ut_date.dt.day,
        "hour_utc": df.ut_date.dt.hour,
        "minute": df.ut_date.dt.minute
    })

    df = pd.concat((date_df, df), axis=1)
    df = df.drop("ut_date", axis=1)

    # Output:
    return loc, df

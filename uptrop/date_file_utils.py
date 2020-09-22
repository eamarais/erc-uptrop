import datetime as dt
import glob
import re
from os import path

import sys


class InvalidSeasonException(Exception):
    pass


class InvalidRegionException(Exception):
    pass


def season_to_date(season):
    if season == "jja":
        start_date = dt.datetime(year=2019, month=6, day=1)
        end_date = dt.datetime(year=2019, month=8, day=31)
    elif season == "son":
        start_date = dt.datetime(year=2019, month=9, day=1)
        end_date = dt.datetime(year=2019, month=11, day=30)
    elif season == "djf":
        start_date = dt.datetime(year=2019, month=12, day=1)
        end_date = dt.datetime(year=2020, month=2, day=29)  # Leap year
    elif season == "mam":
        start_date = dt.datetime(year=2020, month=3, day=1)
        end_date = dt.datetime(year=2020, month=5, day=31)
    elif season == "test":
        start_date = dt.datetime(year=2020, month=3, day=1)
        end_date = dt.datetime(year=2020, month=3, day=3)
    else:
        raise InvalidSeasonException
    return start_date, end_date


def get_tropomi_file_list(trop_dir, date_range):
    """Returns an alphabetically sorted list of Tropomi files
    within a range of dates.

    :param trop_dir: The directory containing the tropomi files
    :type trop_dir: str
    :param date_range: A list of dates. Generation using DateUtil's rrule function is recommended.
    :type date_range: list(datetime)

    :returns: A list of filepaths to tropomi data
    :rtype: list of str
    """
    out = []
    for date in date_range:
        out += (get_tropomi_files_on_day(trop_dir, date))
    return sorted(out)


def get_ocra_file_list(ocra_dir, date_range):
    """Returns an alphabetically sorted list of Ocra files
    within a range of dates.

    :param ocra_dir: The directory containing the ocra files
    :type ocra_dir: str
    :param date_range: A list of dates. Generation using DateUtil's rrule function is recommended.
    :type date_range: list(datetime)

    :returns: A list of filepaths to ocra data
    :rtype: list of str
    """
    out = []
    for date in date_range:
        out += (get_ocra_files_on_day(ocra_dir, date))
    return sorted(out)


def get_tropomi_files_on_day(tomidir, date):
    """Returns a list of tropomi files on a given date.

    Uses the :ref:`get_date` function to extract each candidate file's date from it's filename

    :param tomidir: The directory containing the tropomi files
    :type tomidir: str
    :param date: The date to search for
    :type date: DateTime

    :returns: A list of filepaths to tropomi data
    :rtype: list of str
    """
    # Converts the python date object to a set string representation of time
    # In this case, zero-padded year, month and a datestamp of the Sentinel format
    # See https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
    year = date.strftime(r"%Y")
    month = date.strftime(r"%m")
    datestamp = date.strftime(r"%Y%m%dT")
    tomi_glob_string = path.join(tomidir, 'NO2_OFFL', year, month,'S5P_OFFL_L2__NO2____'+ datestamp + '*')
    tomi_files_on_day = glob.glob(tomi_glob_string)
    print('Found {} tropomi no2 files for {}: '.format(len(tomi_files_on_day), date))
    tomi_files_on_day = sorted(tomi_files_on_day)
    return tomi_files_on_day


def get_ocra_files_on_day(tomidir,date):
    """Returns a list of ocra files on a given date.

    Uses the :ref:`uptrop.tropomi_ut_no2.get_date` function to extract each candidate file's date from it's filename

    :param tomidir: The directory containing the ocra files (usually packaged with tropomi data)
    :type tomidir: str
    :param date: The date to search for
    :type date: DateTime

    :returns: A list of filepaths to ocra data
    :rtype: list of str
    """
    # Get string of day:
    year = date.strftime(r"%Y")
    month = date.strftime(r"%m")
    datestamp = date.strftime(r"%Y%m%dT")
    cld_glob_string = path.join(tomidir, "CLOUD_OFFL", year, month,
                                   'S5P_OFFL_L2__CLOUD__' + datestamp + '*')
    cldfile = glob.glob(cld_glob_string)
    print('Found {} tropomi cloud files for {}: '.format(len(cldfile), date))
    # Order the files:
    cldfile = sorted(cldfile)    
    return cldfile


def get_date(file_name, time_stamp_index = 0):
    """Extracts a datetime object from a filename with a Sentinel timestamp

    See https://regex101.com/r/QNG11l/1 for examples

    :param file_name: The filename to extract the date from
    :type file_name: str
    :param time_stamp_index: Which time-stamp to get the date from if more than one. Defaults to 0.
    :type time_stamp_index: int

    :returns: A DateTime object of the date of the file
    :rtype: DateTime
    """
    # A regular expression that gets Sentinel datestamps out of filenames
    # See https://regex101.com/r/QNG11l/1
    date_regex = r"\d{8}T\d{6}"
    date_string = re.findall(date_regex, file_name)[time_stamp_index]
    # A line for converting Sentinel string reps to datetime
    return dt.datetime.strptime(date_string, r"%Y%m%dT%H%M%S")


def get_gc_files_in_month(gcdir, region, date):
    gc_glob_string = path.join(gcdir, 'nc_sat_files_47L', 'ts_12_15.{}.{}'.format(region, date)+'*')
    # get string of files for this month:
    gcfiles = glob.glob(gc_glob_string)
    
    return sorted(gcfiles)


def get_gc_file_list(gc_dir, region):
    """Gets a list of geoschem files for a given region and set of years

    :param gc_dir: The directory containing the geoschem files
    :type gc_dir: str
    :param region: Can be NA, EU or CH
    :type region: str
    :param date_range: A list of dates. Generation using DateUtil's rrule function is recommended.
    :type date_range: list(datetime)

    :returns: A sorted list of geoschem files
    :rtype: list of str
    """

    # Define target grid:
    if region == 'NA':
        dirreg = '_na'
    elif region == 'EU':
        dirreg = '_eu'
    elif region == 'CH':
        dirreg = '_ch'
    else:
        print("Invalid region; valid regions are 'NA','EU','CH'.")
        raise InvalidRegionException

    # Hard code years and days to process, as it is 2 years of data for
    # one specific season:
    gc_year  = ['2016', '2017']
    gc_month = ['06', '07', '08']

    gc_dir = path.join(gc_dir,'geosfp'+dirreg)
    out_files = []
    for y in gc_year:
        for m in gc_month:
            date = y+m
            day_file = get_gc_files_in_month(gc_dir, region, date)
            # Append to output files:
            for i in range(len(day_file)):
                out_files.append(day_file[i])  
    return sorted(out_files)

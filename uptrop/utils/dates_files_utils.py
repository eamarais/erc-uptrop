import datetime as dt
import calendar
import os
import glob 
import re

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

def season_to_date(season, yyyy):
    if season.lower() == "jja":
        start_month = 6
        end_month = 8
    elif season.lower() == "son":
        start_month = 9
        end_month = 11
    elif season.lower() == "djf":
        start_month = 12
        end_month = 2
    elif season.lower() == "mam":
        start_month = 3
        end_month = 5
    else:
        raise ValueError("Invalid season specified. Accepted seasons are 'djf' or 'DJF', 'mam' or 'MAM', 'jja' or 'JJA', 'son' or 'SON'.")
        
    # Get the start_date
    start_date = dt.datetime(year=yyyy, month=start_month, day=1)
    
    # Get the end_date
    # First check the year: "djf" includes data from this Dec and from Jan+Feb in the next year
    if season == "djf":
        end_yyyy = yyyy+1
    else:
        end_yyyy = yyyy
    
    # Get the last day of the selected yyyy and month using Python's built-in "calendar" module
    last_day_of_end_month = calendar.monthrange(end_yyyy,end_month)[1]
    
    # Get the end_date
    end_date = dt.datetime(year=end_yyyy, month=end_month, day=last_day_of_end_month)
    
    # Get the dates
    return start_date, end_date

def get_files_on_day(data_dir, species, product, date):
    """Returns a list of tropomi files (NO2, O3, HCHO) or cloud files on a given date.

    Uses the :ref:`get_date` function to extract each candidate file's date from it's filename

    :param tomidir: The directory containing the tropomi files
    :type tomidir: str
    :param date: The date to search for
    :type date: DateTime

    :returns: A list of filepaths to data
    :rtype: list of str
    """
    
    # Ensure the species name is valid and use its upper case
    species_upper = species.upper()
    
    if species_upper not in ["NO2", "O3", "HCHO", "CLOUD"]:
        raise ValueError(f"Invalid species: {species}! Available species are 'no2' or 'NO2', 'o3' or 'O3', 'hcho' or 'HCHO', 'cloud' or 'CLOUD'.")
    
    # Ensure the product type is valid
    if product not in ["OFFL", "PAL"]:
        raise ValueError(f"Invalid product type: {product}! Available products types are 'OFFL' and 'PAL'.")
    
    # For "PAL" products, the fine names have one more "_" compared to "OFFL" products
    if product == "PAL":
        prod_file_name = "PAL_"
    else:
        prod_file_name = product
        
    # Create the timestamp for the selected date
    year = date.strftime(r"%Y")
    month = date.strftime(r"%m")
    datestamp = date.strftime(r"%Y%m%dT")
    
    # Decide the "product_name_patterns"
    if species_upper == "NO2": 
        product_name_pattern = f"S5P_{prod_file_name}_L2__{species_upper}____{datestamp}*"
    elif species_upper == "O3":
        product_name_pattern = f"S5P_{prod_file_name}_L2__{species_upper}_____{datestamp}*"
    elif species_upper == "HCHO":
        product_name_pattern = f"S5P_{prod_file_name}_L2__{species_upper}___{datestamp}*"
    elif species_upper == "CLOUD":
        product_name_pattern = f"S5P_{prod_file_name}_L2__{species_upper}__{datestamp}*"
        
    # Search for the files
    glob_string = os.path.join(data_dir, f"{species_upper}_{product}", year, month, product_name_pattern)
    files_on_day = glob.glob(glob_string)
    files_on_day = sorted(files_on_day)
    
    # Notify the user of the files found
    if species == "CLOUD":
        print(f"Found {len(files_on_day)} {species_upper} {product} files for: {date.date()}")
    else:
        print(f"Found {len(files_on_day)} TROPOMI {species_upper} {product} files for: {date.date()}")
    
    # Return the files
    return files_on_day

def get_file_list(data_dir, species, product, date_range):
    """Returns an alphabetically sorted list of files within a range of dates.

    :param ocra_dir: The directory containing the ocra files
    :type ocra_dir: str
    :param date_range: A list of dates. Generation using DateUtil's rrule function is recommended.
    :type date_range: list(datetime)

    :returns: A list of filepaths to ocra data
    :rtype: list of str
    """
    file_list = []
    for date in date_range:
        files_on_the_day = get_files_on_day(data_dir, species, product, date)
        file_list.extend(files_on_the_day)
    return sorted(file_list)
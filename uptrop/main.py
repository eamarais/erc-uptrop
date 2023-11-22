"""
This script provides command-line interface (CLI) tools for the major functions of the 'uptrop' package.

The 'uptrop' package can be used in two ways:

1. As an API (Application Programming Interface) in Python scripts.
For example:
    import uptrop as up
    up.cloud_slice_tropomi(*args, **kwargs)
In this case, the user imports the 'uptrop' module and calls the function `cloud_slice_tropomi` with appropriate arguments.

2. As a CLI tool that can be executed directly from the command line.
For example:
    cloud_slice_tropomi --args1 --args2 --args3
In this case, the user can use the `cloud_slice_tropomi` command directly in the terminal,
passing arguments after the command name.
This method does not require writing a Python script.

Both methods provide the same functionality, and users can choose the one that best suits their workflow.
"""

import os
import logging
import argparse
from .controllers import cloud_slice_tropomi


def parse_args():
    parser = argparse.ArgumentParser(description="Produces a netCDF and preview plot for upper troposphere NO2, O3, or HCHO")
    parser.add_argument("--species", default=None, help="NO2, O3")
    parser.add_argument("--trop_dir", default=None, help="Directory containing TROPOMI data")
    parser.add_argument("--out_dir", default=None, help="Directory to contain finished netCDF4")
    parser.add_argument("--season", default=None, help="Can be jja, son, djf, mam")
    parser.add_argument("--year", default=None, help="Can be 2018, 2019, 2020, 2021, 2022, 2023")
    parser.add_argument("--start_date", default=None, help="Start date of the processing window (yyyy-mm-dd)")
    parser.add_argument("--end_date", default=None, help="End date of the processing window (yyyy-mm-dd)")
    parser.add_argument("--grid_res", default="1x1", help="Can be 05x05, 1x1, 2x25, 4x5")
    parser.add_argument("--cloud_product", default="fresco-wide", help="Can be fresco-wide or o22cld")
    parser.add_argument("--cloud_threshold", default="07", help="Recommended value is 07. Can also test 08, 09, 10")
    parser.add_argument("--fresco_440", default=False, type=bool, help="Use FRESCO-S 440 nm cloud fraction")
    parser.add_argument("--pmin", default=180, type=int, help="Lower bound on cloud height. Defaults to 180.")
    parser.add_argument("--pmax", default=450, type=int, help="Upper bound on cloud height. Defaults to 450.")
    parser.add_argument("--product", default=None, help="TROPOMI product name. Can be OFFL or PAL")
    parser.add_argument("--version", default="v1", help="Version number to append to filename")
    parser.add_argument("--save_data", default=True, type=bool, help="Save data or not")
    args = parser.parse_args()
    return args


def main_tropomi():
    # Use the command line arguments to configure your program
    args = parse_args()

    # Define log file name based on args
    log_file_name = (
    f'tropomi_{args.species}_{args.start_date}_{args.end_date}_'
    f'{args.grid_res}_{args.cloud_product}_{args.cloud_threshold}_{args.pmin}_'
    f'{args.pmax}_{args.product}_{args.version}.log'
    )

    log_file = os.path.join(args.out_dir, log_file_name)

    # Set up the logger
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M')

    cloud_slice_tropomi(
        species=args.species,
        trop_dir=args.trop_dir,
        out_dir=args.out_dir,
        season=args.season,
        year=args.year,
        start_date=args.start_date,
        end_date=args.end_date,
        grid_res=args.grid_res,
        cloud_product=args.cloud_product,
        cloud_threshold=args.cloud_threshold,
        fresco_440=args.fresco_440,
        pmin=args.pmin,
        pmax=args.pmax,
        product=args.product,
        version=args.version,
    )

# Code below here will only be excuted when this script is used directly, not as a module
# This allows you to test the CLI tools while developing the software
if __name__ == "__main__":
    main_tropomi()


"""
Example usages:
export PYTHONPATH=$PYTHONPATH:/lustre/scratch/scratch/ucfauxx/projects/uptrop/erc-uptrop/

python -m uptrop.main --species "no2" --product "PAL" --trop_dir "/lustre/projects/uptrop/tropomi/Data/" --out_dir "/lustre/scratch/scratch/ucfauxx/projects/uptrop/scratch/try-v2" --start_date "2020-08-31" --end_date "2020-08-31" --pmin 450 --pmax 600
python -m uptrop.main --species "o3" --product "OFFL" --trop_dir "/lustre/projects/uptrop/tropomi/Data/" --out_dir "/lustre/scratch/scratch/ucfauxx/projects/uptrop/scratch/try-v2" --start_date "2020-08-31" --end_date "2020-08-31" --pmin 450 --pmax 600

"""
import sys, os
# Import hack
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..'))

from uptrop import date_file_utils
import datetime as dt
from dateutil import rrule as rr

TROPOMI_DIR="/data/uptrop/nobackup/tropomi/Data/"
PANDORA_DIR="/data/uptrop/nobackup/pandora/"
GC_DIR="/data/uptrop/Projects/DEFRA-NH3/GC/"
START_DATE=dt.datetime(year=2016, month=6, day=1)
END_DATE=dt.datetime(year=2016, month=6, day=3)

def test_get_gc_file_list():
    date_range = rr.rrule(rr.DAILY, dtstart=START_DATE, until=END_DATE)
    out = date_file_utils.get_gc_file_list(GC_DIR, 'NA', date_range)
    assert len(out) > 1

#!/usr/bin/python

import os
import sys
import numpy as np

# Import hack
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..'))
from uptrop.constants import G as g
from uptrop.constants import MW_AIR as mmair
from uptrop.constants import R_GAS_UNIV as rgas

def alt2pres(height):
    """Convert altitude provided in m to pressure in Pa.

    :param height: Height in m.

    :return: Pressure in Pa.
    """

    # Define additional constants:
    p_zero=101325    #Pa
    lapse=6.49       #K/km
    t_zero=288.15    #K

    # Convert lapse rate frrom K/km to K/m:
    lapse=np.multiply(lapse,1e-3)

    # Calculate pressure:
    tval1=np.divide((np.multiply(g,mmair)),np.multiply(rgas,lapse))
    tval2=np.divide(np.multiply(lapse,height),t_zero)
    pressure=np.multiply(p_zero,(np.power(np.subtract(1.0,tval2),tval1)))    

    # Output goes here:
    return (pressure)

 

import numpy as np
from .constants import G as g
from .constants import MW_AIR as mmair
from .constants import R_GAS_UNIV as rgas

import logging

def alt2pres(height):
    """ 
    When the height is big enough, 1-tval2 < 0, then this raises
    RuntimeWarning: invalid value encountered in power
      pressure=np.multiply(p_zero,(np.power(np.subtract(1.0,tval2),tval1)))

    """
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
    logging.info(f"check 1: {np.subtract(1.0,tval2)}")
    logging.info(f"check 2: {tval1}")
    pressure=np.multiply(p_zero,(np.power(np.subtract(1.0,tval2),tval1)))    

    # Output:
    logging.info(f"pressure: {pressure}")
    return (pressure)

def pres2alt(pressure):
    """Convert pressure provided in Pa to altitude in m.

    :param pressure: Pressure in Pa.

    :return: Height in m.
    """

    # Define additional constants:
    p_zero=101325    #Pa

    # Calculate height:
    height=(1 - ((pressure/p_zero)**0.190284)) * 145366.45 * 0.3048

    # Output:
    return (height)
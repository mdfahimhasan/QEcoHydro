# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

import numpy as np


def calc_esat(Tc):
    """
    Calculate saturation vapor pressure from temperature (Eq. 2.17 in TH)
    :param Tc: Temperature (deg C).
    :return: Saturation vapor pressure (kPa).
    """
    esat = 0.6108 * np.exp(17.27 * Tc / (Tc + 237.3))  # Eq. 2.17 TH

    return esat


def calc_VPD(esat, e):
    """
    Calculate vapor pressure deficit (Eq. 2.10 in TH)
    :param esat: Saturation vapor pressure (kPa).
    :param e: Vapor pressure (kPa).

    :return: Vapor pressure deficit (kPa).
    """
    VPD = esat - e  # Eq 2.10 TH

    return VPD


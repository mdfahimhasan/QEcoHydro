# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

import numpy as np


def calc_lambda(Temp):
    """
    Calculate latent heat of vaporization (lambda) at given temperature.

    :param Temp: Temperature (deg celsius).

    :return: Latent heat of vaporization (kJ/kg).
    """
    # Latent Heat of Vaporization (kJ/kg)
    lambda_ = 1000 * (2.501 - 0.002361 * Temp)  # Latent Heat of Vaporization (kJ/kg); TH eq. 2.1

    return lambda_


def calc_Lnet(SW_in, Rso, vp, temp, region='arid'):
    """
    Calculate net longwave radiation.

    :param SW_in: Incoming solar radiation (W/m2).
    :param Rso: Clear sky solar radiation (W/m2).
    :param vp: Vapor pressure (kpa).
    :param temp: Temperature (deg C).
    :param region: either 'arid' or 'humid'. Affects calculation of f.

    :returns: net longwave radiation (W/m2).
    """
    # empirical cloud factor
    # SW_in is calculated by TH eq. 5.16 or given.
    # Rso is calculated using dates, lat, elev (look into Ra_Rso.py for detail)
    if region == 'arid':
        f = 1.35 * (SW_in / Rso) - 0.35
    else:  # for humid
        f = SW_in / Rso  # unitless

    # effective emissivity
    e_prime = 0.34 - 0.14 * np.sqrt(vp / 1000)  # TH eq. 5.23; unitless

    # net longwave radiation
    stf_boltzman_const = 5.67 * 10 ** -8  # unit W m-2 K-4
    L_net = -f * e_prime * stf_boltzman_const * (temp + 273.15) ** 4  # TH eq. 5.22

    return L_net


def calc_Rn_RC(SW_in, Lnet):
    """
    Calculate net radiation for reference crop.

    :param SW_in: Incoming shortwave radiation (W/m2).
    :param Lnet: Net longwave radiation (W/m2).
    :param temp: Temperature (deg C).

    :return: Net solar radiation in W m-2.
    """
    # Net Radiation over a Reference Crop, in W/m2
    # albedo considered 0.23
    albedo = 0.23
    Rn_RC = SW_in * (1 - albedo) + Lnet  # TH eq. 5.28

    return Rn_RC

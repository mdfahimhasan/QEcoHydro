# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

import numpy as np
from vars_constants import EnvConstants, ForestConstants, GrassConstants


# # Aerodynamic Resistance Parameterization
def calc_ra(wind_spd, canopy_type='grass', params=None):
    """
    Calculates aerodynamic resistance around a canopy.

    :param wind_spd: Wind speed (m/s).
    :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the
                          vars_constants.py file. Default set to 'grass'.
    :param params: A dictionary of params zm, d, z0 to adjust manually.
                 Default set to None to use these vars from the vars_constants.py.

    :return: Aerodynamic resistance (s m-1).
    """
    global zm, d, z0

    k = EnvConstants.k  # von Karman's constant

    if canopy_type == 'grass':
        zm = GrassConstants.zm  # Measurement height for wind speed (m)
        d = GrassConstants.d  # Zero plane displacement height (m)
        z0 = GrassConstants.z0  # Aerodynamic roughness height (m)

    elif canopy_type == 'forest':
        zm = ForestConstants.zm  # Measurement height for wind speed (m)
        d = ForestConstants.d  # Zero plane displacement height (m)
        z0 = ForestConstants.z0  # Aerodynamic roughness height (m)

    # Using user defined parameters if params is not None
    if params is not None:
        zm = params['zm']
        d = params['d']
        z0 = params['z0']

    # assuming that measurement of height of temperature and wind speed are the same
    ra = 1 / (k ** 2 * wind_spd) * np.log((zm - d) / z0) * np.log((zm - d) / (z0 / 10))  # TH Eq 22.9

    return ra


# # Surface Resistance Parameterization
# Using the Jarvis-Stewart scheme (ref: Terrestrial Hydrometry by Shuttleworth)

def calc_rs(SW_in, vpd, temp, FC, WP, Su_0, canopy_type='grass', params=None):
    """
    Calculates stomatal resistance.

    :param SW_in: Incoming shortwave radiation (W m-2).
    :param vpd: Vapor pressure deficit (kPa).
    :param temp: Temperature (in deg C).
    :param FC: Field capacity.
    :param WP: Wilting point.
    :param Su_0: Updated unsaturated zone storage (Soil moisture (mm)) after subtracting recharge.
    :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the
                          vars_constants.py file. Default set to 'grass'.
    :param params: A dictionary of params g0 and gc to adjust manually.
         Default set to None to use these vars from the vars_constants.py.

    :return: Total surface resistance (s m-1).
    """
    global g0, gc

    # # atmospheric and  canopy constants
    # for gR
    KR = EnvConstants.KR

    # for gD
    KD1 = EnvConstants.KD1
    KD2 = EnvConstants.KD2

    # for gT
    TL = EnvConstants.TL
    T0 = EnvConstants.T0
    TH = EnvConstants.TH
    aT = EnvConstants.aT

    # for gSM, gS
    if canopy_type == 'grass':
        g0 = GrassConstants.g0
        gc = EnvConstants.gc

    elif canopy_type == 'forest':
        g0 = ForestConstants.g0
        gc = EnvConstants.gc

    # Using user defined parameters if params is not None
    if params is not None:
        g0 = params['g0']
        gc = params['gc']

    def calc_gR(SW_in):
        """
        Calculate resistance to soil heat flux (Eq. 24.2 in TH)
        :param SW_in: Incoming shortwave radiation (W m-2).
        :return: Resistance to soil heat flux (none).
        """
        gR = SW_in * (1000 + KR) / (1000 * (SW_in + KR))  # Eq. 24.2 TH
        gR = max(gR, 0.)

        return gR

    def calc_gD(VPD):
        """
        Calculate resistance to vapor pressure deficit (Eq. 24.3 in TH)
        :param VPD: Vapor pressure deficit (kPa).
        :return: Resistance to vapor pressure deficit (none).
        """
        gD = 1 + KD1 * VPD + KD2 * VPD ** 2  # Eq. 24.3 TH
        gD = max(gD, 0.)

        return gD

    def calc_gT(Tc):
        """
        Calculate resistance to temperature (Eq. 24.4 in TH)
        :param Tc: Temperature (deg C).

        :return: Resistance to temperature (none).
        """
        Tk = Tc + 273.17  # unit from deg C to Kelvin
        gT = ((Tk - TL) * ((TH - Tk) ** aT)) / ((T0 - TL) * ((TH - T0) ** aT))  # Eq. 24.4 TH, problem statement B.1

        gT = max(gT, 0.)

        return gT

    def calc_gSM(FC, WP, Su_0):
        """
        Calculate resistance to soil moisture flux (Eq. 24.6 in TH)
        :param FC: Field capacity.
        :param WP: Wilting point.
        :param Su_0: Updated unsaturated zone storage (Soil moisture (mm)) after subtracting recharge.

        :return: Resistance to soil moisture flux (none).
        """
        gSM = min(max((Su_0 - WP) / (FC - WP), 0), 1)

        return gSM

    def calc_gS(g0, gc, gR, gD, gT, gSM):
        """
        Calculate total surface conductance (Eq. 24.1 in TH)

        :return: Total surface conductance (mm s-1).
        """
        gS = g0 * gc * gR * gD * gT * gSM  # Eq. 24.1 TH
        gS = max(gS, 0.0001)

        return gS

    # calculating stress factors (resistance)
    gR = calc_gR(SW_in=SW_in)
    gD = calc_gD(VPD=vpd)
    gT = calc_gT(Tc=temp)
    gSM = calc_gSM(FC=FC, WP=WP, Su_0=Su_0)

    # calculating gS
    gS = calc_gS(g0=g0, gc=gc, gR=gR, gD=gD, gT=gT, gSM=gSM)

    # calculating stomatal resistance
    rs = 1 / gS  # Eq. 24.1 TH, unit s mm-1
    rs = rs * 1000  # unit in s m-1

    rs = min(rs, 10 ** 6)

    return rs

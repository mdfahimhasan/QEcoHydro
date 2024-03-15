import numpy as np
from vars_constants import EnvConstants


class atmosphere:
    """
    Calculates atmospheric variables like vapor pressure, saturated vapour pressure, and vapor pressure deficit.

    :param P: Atmospheric pressure (kPa).

    :param e: Vapor pressure (kPa). Defaults set to None.
    :param esat: Saturation vapor pressure (kPa). Defaults set to None.
    :param VPD: Vapor pressure deficit (kPa).  Defaults set to None.

    """

    def __init__(self):
        # atmospheric and  canopy constants
        self.P = EnvConstants.P

        # variables estimated using methods (functions) within the class
        self.e = None
        self.esat = None
        self.VPD = None

    def calc_e(self, q):
        """
        Calculate vapor pressure from specific humidity and pressure (Eq. 2.9 in TH)
        :param q: Specific humidity (g kg-1).
        :param P: Atmospheric pressure (kPa) from EnvConstants from "vars.py".
        :return: Vapor pressure (kPa).
        """
        q = q / 1000  # converting g/kg to kg/kg

        self.e = self.P * q / 0.622  # Eq 2.9 TH
        return self.e

    def calc_esat(self, Tc):
        """
        Calculate saturation vapor pressure from temperature (Eq. 2.17 in TH)
        :param Tc: Temperature (deg C).
        :return: Saturation vapor pressure (kPa).
        """
        self.esat = 0.6108 * np.exp(17.27 * Tc / (Tc + 237.3))  # Eq. 2.17 TH

        return self.esat

    def calc_VPD(self, esat, e):
        """
        Calculate vapor pressure deficit (Eq. 2.10 in TH)
        :param esat: Saturation vapor pressure (kPa).
        :param e: Vapor pressure (kPa).

        :return: Vapor pressure deficit (kPa).
        """
        self.VPD = esat - e  # Eq 2.10 TH

        return self.VPD


# This block will help to check the functions if the .py is run as a script (not imported as module)
if __name__ == '__main__':
    print('Calculating vapor pressure deficit value')
    atm = atmosphere()
    e = atm.calc_e(2.444)
    print(f'e = {e} Kpa')
    e_sat = atm.calc_esat(12.697)
    print(f'e_sat = {e_sat} Kpa')
    D = atm.calc_VPD(esat=e_sat, e=e)
    print(f'VPD = {D} Kpa')
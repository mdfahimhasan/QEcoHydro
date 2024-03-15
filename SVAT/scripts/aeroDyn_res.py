import numpy as np
from vars_constants import EnvConstants, ForestConstants, GrassConstants


# Part A: Aerodynamic Parameterization

class aerodynamic_res:

    def __init__(self, canopy_type='grass'):
        """
        Calculates aerodynamic resistance around a canopy.

        :param grass_or_forest: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the
                                vars.py file. Default set to 'grass'.

        :param k (float): von Karman's constant (none). Comes from "vars.py".
        :param h (float): Canopy height (m). Comes from "vars.py" based on 'grass'/'forest'.
        :param LAI (float): Leaf area index (none). Comes from "vars.py" based on 'grass'/'forest'.
        :param zm (float): Measurement height for wind speed (m). Comes from "vars.py" based on 'grass'/'forest'.
        :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the
                                vars.py file. Default set to 'grass'.

        :param d: Zero plane displacement height (m). Defaults set to None.
        :param z0: Aerodynamic roughness of crop (m). Defaults set to None.
        :param ra: Aerodynamic resistance (s m-1). Defaults set to None.
        """
        # atmospheric and  canopy constants
        self.k = EnvConstants.k

        if canopy_type == 'grass':
            self.h = GrassConstants.h
            self.LAI = GrassConstants.LAI
            self.zm = GrassConstants.zm
        elif canopy_type == 'forest':
            self.h = ForestConstants.h
            self.LAI = ForestConstants.LAI
            self.zm = ForestConstants.zm

        # variables estimated using methods (functions) within the class
        self.d = None
        self.z0 = None
        self.ra = None

    def calc_d(self):
        """
        Calculate zero plane displacement height using constants from either ForestConstants or GrassConstants. (TH Eq 22.2)
        :return: Zero plane displacement height (m).
        """
        self.d = 1.1 * self.h * np.log(1 + (self.LAI / 5) ** 0.25)  # TH Eq 22.2

        return self.d

    def calc_z0(self):
        """
        Calculate aerodynamic roughness of crop. (TH Eq 22.3, 22.4)
        :return: Aerodynamic roughness of crop (m).
        """
        self.z0 = 0.3 * self.h * (1 - self.d / self.h)  # TH Eq 22.4, for LAI > 1

        return self.z0

    def calc_ra(self, wind_spd):
        """
        Calculate aerodynamic resistance using constants from Environmental crop_constants. (TH Eq 22.9)
        :param wind_spd: Wind speed (m/s).

        :return: Aerodynamic resistance (s m-1).
        """
        # assuming that measurement of height of temperature and wind speed are the same
        self.ra = 1 / (self.k ** 2 * wind_spd) * np.log((self.zm - self.d) / self.z0) * np.log((self.zm - self.d) / (self.z0 / 10))  # TH Eq 22.9

        return self.ra


# This block will help to check the d, z0 values if the .py is run as a script (not imported as module)
if __name__ == '__main__':
    # Checking values for grass
    print('=== Grass ===')
    aero_dyn_g = aerodynamic_res(canopy_type='grass')
    print('d =', aero_dyn_g.calc_d())
    print('z0 =', aero_dyn_g.calc_z0())
    print('ra =', aero_dyn_g.calc_ra(wind_spd=2.43), '(s m-1)')

    # Checking values for forest
    print('=== Forest ===')
    aero_dyn_f = aerodynamic_res(canopy_type='forest')
    print('d =', aero_dyn_f.calc_d())
    print('z0 =', aero_dyn_f.calc_z0())
    print('ra =', aero_dyn_f.calc_ra(wind_spd=2.43), '(s m-1)')

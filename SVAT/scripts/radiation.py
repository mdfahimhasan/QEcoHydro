# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

from vars_constants import EnvConstants, ForestConstants, GrassConstants


# Part C: Radiation Parameterization

class radiation:

    def __init__(self, canopy_type='grass'):
        """
        Calculates upward longwave radiation and net radiation.

        :param cp: Specific heat capacity of air (J kg-1 K-1). Environmental Constant.
        :param rho_a: Moist Air density (kg m-3). Environmental Constant.
        :param sigma: Stefan-Boltzmann constant (W m-2 K-4).
        :param emissivity: Emissivity of grass/forest (none).
        :param a: Albedo (none). Crop Constant.
        :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the
                            vars.py file. Default set to 'grass'.

        :param LW_up: Long wave radiation from surface (W m-2). Defaults set to None.
        :param Rn: Net solar radiation (W m-2). Defaults set to None.
        """
        # atmospheric and  canopy constants
        self.cp = EnvConstants.cp
        self.rho_a = EnvConstants.rho_a
        self.sigma = EnvConstants.sigma

        if canopy_type == 'grass':
            self.emissivity = GrassConstants.E_surface
            self.a = GrassConstants.a
            self.LW_up_init = GrassConstants.LW_up_init

        elif canopy_type == 'forest':
            self.emissivity = ForestConstants.E_surface
            self.a = ForestConstants.a
            self.LW_up_init = ForestConstants.LW_up_init

        # variables estimated using methods (functions) within the class
        self.LW_up = None
        self.Rn = None

    def calc_LW_up(self, Ts1_k, Ts2_k):
        """
        Calculate long wave radiation from surface temperature and emissivity (Eq 5.19 in TH).

        ** This function will start to be implemented from the 3rd model step because initially we don't know values of
        Ts (as we don't know H). In the 1st and 2nd steps, values are chosen from typical values calculated for the next
        few hours. Once we calculate Ts following 'Canopy Water Balance Parameterization' module, this function can be
        implemented to calculate LW_up.

        :param Ts1_k , Ts2_k: Surface temperature at previous 2 steps (Kelvin).
        :param emissivity: Emissivity of crop (none). Environmental constant.
        :param sigma: Stefan-Boltzmann constant (W m-2 K-4). Environmental constant.

        :return: Longwave radiation from surface (W m-2).
        """

        LW_up = - (self.emissivity * self.sigma * (Ts1_k ** 4 + Ts2_k ** 4) / 2)  # Modified Eq 5.19 TH (Problem statement C.2)
                                                                                # (-ve) sign for upward flux
        self.LW_up = LW_up

        return self.LW_up

    def calc_Rn(self, SW_in, LW_up, LW_down):
        """
        Calculate net radiation (Eq 5.28 in TH).
        :param SW_in: Incoming shortwave radiation (W m-2).
        :param LW_up: Long wave radiation from surface (W m-2).
        :param LW_down: Long wave radiation from atmosphere (W m-2).

        :return: Net solar radiation (W m-2).
        """
        LW_net = LW_down + LW_up  # LW_up will always be (-ve) as upward energy flux
        Rn = (1 - self.a) * SW_in + LW_net  # Eq 5.28 TH
        self.Rn = Rn

        return self.Rn


# This block will help to check the functions if the .py is run as a script (not imported as module)
if __name__ == '__main__':
    # Checking values for grass
    print('=== Grass ===')
    rd_g = radiation(canopy_type='grass')
    print(f'SWin = {500}')
    lw_up = rd_g.calc_LW_up(Ts1_k=293, Ts2_k=293.5)
    print('LW up = ', lw_up)
    print('Rn = ', rd_g.calc_Rn(SW_in=500, LW_up=lw_up, LW_down=265))

    # # Checking values for forest
    print('=== Forest ===')
    rd_f = radiation(canopy_type='forest')
    print(f'SWin = {500}')
    lw_up = rd_f.calc_LW_up(Ts1_k=293, Ts2_k=293.5)
    print('LW up = ', lw_up)
    print('Rn = ', rd_f.calc_Rn(SW_in=500, LW_up=lw_up, LW_down=265))

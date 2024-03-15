import numpy as np
from vars_constants import EnvConstants, ForestConstants, GrassConstants


# Part B: Surface Resistance Parameterization

class surface_res:

    def __init__(self, canopy_type='grass'):
        """
        Calculates resistance to soil heat flux (gR), vapor pressure deficit (gD), temperature (gT),
        soil moisture flux (gSM), total surface conductance (gS), and
        the final surface resistance (rs) in a canopy.

        :param KR: Parameter in Equ (24.2) (W m-2). Environmental constant.
        :param KD1: Parameter in Equ (24.6) (kPa-1). Environmental constant.
        :param KD2: Parameter in Equ (24.6) (kPa-2). Environmental constant.
        :param TL: Parameter in Equ (24.4) and Equ (24.5) (K). Environmental constant.
        :param T0: Parameter in Equ (24.4) and Equ (24.5) (K). Environmental constant.
        :param TH: Parameter in Equ (24.4) and Equ (24.5) (K). Environmental constant.
        :param aT: Parameter in temperature stress factor (none). Environmental constant.
        :param KM2: Parameter in Equ (24.6) (mm-1). Environmental constant.
        :param g0: Canopy Specific Constant (mm s-1).
        :param gc: Canopy cover factor (none).
        :param KM1: Parameter in Equ (24.6) (none). Environmental constant.
        :param SMo: Maximum soil moisture accessible to roots. Parameter in Equ (24.6) (mm).
        :param SMinit: Initial Root-accessible soil moisture (mm).
        :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the
                             vars.py file. Default set to 'grass'.

        :param gR: Resistance to soil heat flux (none). Defaults set to None.
        :param gD: Resistance to vapor pressure deficit (none). Defaults set to None.
        :param gT: Resistance to temperature (none). Defaults set to None.
        :param gSM: Resistance to soil moisture flux (none). Defaults set to None.
        :param gS: Total surface conductance (mm s-1). Defaults set to None.
        :param rs: Total surface resistance (s m-1). Defaults set to None.
        """
        # # atmospheric and  canopy constants
        # for gR
        self.KR = EnvConstants.KR

        # for gD
        self.KD1 = EnvConstants.KD1
        self.KD2 = EnvConstants.KD2

        # for gT
        self.TL = EnvConstants.TL
        self.T0 = EnvConstants.T0
        self.TH = EnvConstants.TH
        self.aT = EnvConstants.aT

        # for gSM
        self.KM2 = EnvConstants.KM2

        # for gSM, gS
        if canopy_type == 'grass':
            self.g0 = GrassConstants.g0
            self.gc = EnvConstants.gc
            self.KM1 = GrassConstants.KM1
            self.SMo = GrassConstants.SMo
            self.SMinit = GrassConstants.SMinit

        elif canopy_type == 'forest':
            self.g0 = ForestConstants.g0
            self.gc = EnvConstants.gc
            self.KM1 = ForestConstants.KM1
            self.SMo = ForestConstants.SMo
            self.SMinit = ForestConstants.SMinit

        # variables estimated using methods (functions) within the class
        self.gR = None
        self.gD = None
        self.gT = None
        self.gSM = None
        self.gS = None

        self.rs = None

    def calc_gR(self, SW_in):
        """
        Calculate resistance to soil heat flux (Eq. 24.2 in TH)
        :param SW_in: Incoming shortwave radiation (W m-2).
        :return: Resistance to soil heat flux (none).
        """
        gR = SW_in * (1000 + self.KR) / (1000 * (SW_in + self.KR))  # Eq. 24.2 TH
        self.gR = max(gR, 0.)

        return self.gR

    def calc_gD(self, VPD):
        """
        Calculate resistance to vapor pressure deficit (Eq. 24.3 in TH)
        :param VPD: Vapor pressure deficit (kPa).
        :return: Resistance to vapor pressure deficit (none).
        """
        gD = 1 + self.KD1 * VPD + self.KD2 * VPD ** 2  # Eq. 24.3 TH
        self.gD = max(gD, 0.)

        return self.gD

    def calc_gT(self, Tc):
        """
        Calculate resistance to temperature (Eq. 24.4 in TH)
        :param Tc: Temperature (deg C).

        :return: Resistance to temperature (none).
        """
        Tk = Tc + 273.17  # unit from deg C to Kelvin
        gT = ((Tk - self.TL) * ((self.TH - Tk) ** self.aT)) / (
                    (self.T0 - self.TL) * ((self.TH - self.T0) ** self.aT))  # Eq. 24.4 TH, problem statement B.1

        self.gT = max(gT, 0.)

        return self.gT

    def calc_gSM(self, SM):
        """
        Calculate resistance to soil moisture flux (Eq. 24.6 in TH)
        :param SM: Soil moisture (mm). At the beginning of time step, this value have to be set to SMinit.
                    After that, at each time-step, SM will be set as the SMnew of the previous step.

        :return: Resistance to soil moisture flux (none).
        """
        gSM = 1 - self.KM1 * np.exp(self.KM2 * (SM - self.SMo))  # Eq. 24.6 TH
        self.gSM = max(gSM, 0.)

        return self.gSM

    def calc_gS(self):
        """
        Calculate total surface conductance (Eq. 24.1 in TH)
        :return: Total surface conductance (mm s-1).
        """
        gS = self.g0 * self.gc * self.gR * self.gD * self.gT * self.gSM  # Eq. 24.1 TH
        self.gS = max(gS, 0.0001)

        return self.gS

    def calc_rs(self, gS):
        """
        Calculate total surface resistance (Eq. 24.1 in TH)
        :return: Total surface resistance (s m-1).
        """
        rs = 1 / gS  # Eq. 24.1 TH, unit s mm-1
        rs = rs * 1000  # unit in s m-1

        self.rs = min(rs, 10 ** 6)

        return self.rs


# This block will help to check the functions if the .py is run as a script (not imported as module)
if __name__ == '__main__':
    # Checking values for grass
    print('=== Grass ===')
    sr_g = surface_res(canopy_type='grass')
    print(f'SWin = {100}')
    print('gR = ', sr_g.calc_gR(SW_in=99))
    print('gD = ', sr_g.calc_gD(VPD=1.070))
    print('gT = ', sr_g.calc_gT(Tc=10.68))
    print('gSM = ', sr_g.calc_gSM(SM=10))
    print('g0 = ', sr_g.g0)
    print('gS = ', sr_g.calc_gS())
    print('rS = ', sr_g.calc_rs(gS=sr_g.gS), '(s m-1)')

    # Checking values for forest
    print('=== Forest ===')
    sr_f = surface_res(canopy_type='forest')
    print(f'SWin = {100}')
    print('gR = ', sr_f.calc_gR(SW_in=99))
    print('gD = ', sr_f.calc_gD(VPD=1.070))
    print('gT = ', sr_f.calc_gT(Tc=10.68))
    print('gSM = ', sr_f.calc_gSM(SM=10))
    print('g0 = ', sr_f.g0)
    print('gS = ', sr_f.calc_gS())
    print('rS = ', sr_f.calc_rs(gS=sr_g.gS), '(s m-1)')


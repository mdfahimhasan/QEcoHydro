# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

from vars_constants import EnvConstants, ForestConstants, GrassConstants


# Part D: Canopy Water Balance Parameterization (Includes Interception)

class canopy_WB:

    def __init__(self, canopy_type='grass'):
        """
        Calculates Latent heat, delta, psychrometric constant, canopy storage fluxes, canopy drainage,
        evaporation, transpiration, total ET, sensible heat, and surface temperature at each step.

        :param cp: specific heat of dry air at constant pressure (kJ/kg/K).
        :param P: pressure (kPa).
        :param rho_a: Moist Air density (kg m-3).
        :param S: Maximum Canopy Water Storage (mm).
        :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the
                                 vars.py file. Default set to 'grass'.

        :param LH: Latent heat of vaporization (lambda) (J/kg). Defaults set to None.
        :param delta: Slope of the saturation vapor pressure curve (kPa/°C). Defaults set to None.
        :param psy_const: Psychrometric constant (kPa/°C). Defaults set to None.
        :param Cinterim: Interim canopy storage (mm). Defaults set to None.
        :param Cactual: Actual canopy storage (mm). Defaults set to None.
        :param Dcanopy: Canopy drainage (mm). Defaults set to None.
        :param lambdaEi: Evaporation rate of intercepted rainfall with free water canopy (W m-2). Defaults set to None.
        :param Cfinal: Final canopy storage (mm). Defaults set to None.
        :param lambdaEt: Transpiration rate from canopy (W m-2). Defaults set to None.
        :param lambdaET: Total evapotranspiration (ET) rate (W m-2). Defaults set to None.
        :param H: Sensible heat flux (J/kg). Defaults set to None.
        :param Ts_k: Surface temperature (K). Defaults set to None.
        """
        # atmospheric and  canopy constants
        self.cp = EnvConstants.cp
        self.P = EnvConstants.P
        self.rho_a = EnvConstants.rho_a

        if canopy_type == 'grass':
            self.S = GrassConstants.S

        elif canopy_type == 'forest':
            self.S = ForestConstants.S

        # variables estimated using methods (functions) within the class
        self.LH = None
        self.delta = None
        self.psy_const = None

        self.Cinterim = None
        self.Cactual = None
        self.Dcanopy = None
        self.lambdaEi = None
        self.Cfinal = None
        self.lambdaEt = None
        self.lambdaET = None
        self.H = None
        self.Ts_k = None

    def calc_LH(self, Tc):
        """
        Calculate latent heat of vaporization from temperature (Eq 2.1 in TH)
        :param Tc: temperature (in deg C).

        :return: latent heat of vaporization (lambda) (J/kg).
        """
        self.LH = (2.501 - 0.002361 * Tc) * 10 ** 6  # Eq 2.1 TH, problem statement D.1.1

        return self.LH

    def calc_delta(self, Tc, e_sat):
        """
        Calculate gradient of the saturation vapor pressure curve from saturation vapor pressure and temperature
        (Eq 2.18 in TH).
        :param e_sat: saturation vapor pressure  (kPa).
        :param Tc: temperature in (deg C).

        :return: slope of the saturation vapor pressure curve (kPa/°C).
        """
        self.delta = 4098 * e_sat / (Tc + 237.3) ** 2  # Eq 2.18 TH, problem statement D.1.1

        return self.delta

    def calc_psy_const(self, LH):
        """
        Calculate psychrometric constant (Eq 2.25 in TH).
        :param LH: latent heat of vaporization (lambda) (J/kg).

        :return: psychrometric constant (kPa/°C).
        """
        self.psy_const = self.cp * self.P / (0.622 * LH)  # Eq. 2.25 TH, problem statement D.1.1

        return self.psy_const

    def calc_Cinterim(self, precip, Cfinal_prev):
        """
        Calculate interim canopy drainage (D.1.6 in problem statement).
        For first time step, it has to be set to zero.
        :param precip: precipitation (mm).
        :param Cfinal_prev: final canopy storage (mm) from previous step.

        :return: interim canopy storage (mm).
        """
        Cinterim = precip + Cfinal_prev  # Given in problem statement D.1.6
        self.Cinterim = max(Cinterim, 0)  # Ensuring positive value

        return self.Cinterim

    def calc_Cactual(self):
        """
        Calculate actual canopy storage for each row.

        :return: actual canopy storage (mm).
        """
        if self.Cinterim > self.S:
            self.Cactual = self.S
        else:
            self.Cactual = self.Cinterim

        return self.Cactual

    def calc_Dcanopy(self):
        """
        Calculate canopy drainage for each row.

        :return: canopy drainage (mm).
        """
        if self.Cinterim > self.S:
            self.Dcanopy = self.Cinterim - self.S
        else:
            self.Dcanopy = 0

        return self.Dcanopy

    def calc_lambdaEi(self, ra, Rn, VPD_ref, Cactual):
        """
        Calculate evaporation rate of intercepted rainfall with free water canopy (Eq 22.14 in TH).

        :param ra: aerodynamic resistance (s/m).
        :param Rn: net radiation (W/m^2).
        :param VPD_ref: vapor pressure deficit (kPa).
        :param Cactual: actual canopy storage (mm).

        :return: evaporation rate of intercepted rainfall with free water canopy (W m-2).
        """
        A = Rn  # Available energy = net radiation problem statement D.1.8

        self.lambdaEi = (Cactual / self.S) * (self.delta * A + (self.rho_a * self.cp * VPD_ref / ra)) / \
                        (self.delta + self.psy_const)  # Eq 22.14 TH, problem statement D.1.4

        return self.lambdaEi

    def calc_Cfinal(self):
        """
        Calculate final canopy storage.

        :return: final canopy storage (mm).
        """
        m_to_mm = 1000
        rho_w = 1000  # Density of water in kg/m^3

        lambdaEi_mm = (self.lambdaEi * 3600 / self.LH / rho_w) * m_to_mm  # Convert from W m-2 to mm

        self.Cfinal = max(self.Cactual - lambdaEi_mm, 0)  # Problem statement D.1.9

        return self.Cfinal

    def calc_lambdaEt(self, ra, rs, Rn, VPD_ref):
        """
        Calculate transpiration rate (Eq 22.18 in TH).

        :param ra: aerodynamic resistance (s/m).
        :param rs: canopy surface resistance (s/m).
        :param Rn: net radiation (W/m^2).
        :param VPD_ref: vapor pressure deficit (kPa).

        :return: transpiration rate from canopy (W m-2).
        """
        A = Rn  # Available energy = net radiation problem statement D.1.8

        self.lambdaEt = (self.delta * A + (self.rho_a * self.cp * VPD_ref / ra)) / \
                        (self.delta + self.psy_const * (1 + (rs/ra)))  # Eq 22.14 TH, problem statement D.1.4

        return self.lambdaEt

    def calc_lambdaET(self):
        """
        Calculate total evaporation rate (Eq 22.17 in TH).

        :return: total evapotranspiration (ET) rate (W m-2).
        """
        # the 2nd part of the equation with the (1-C/s) term is the portion coming through (from soil moisture) that
        # is used for canopy transpiration
        self.lambdaET = self.lambdaEi + (1 - self.Cactual / self.S) * self.lambdaEt  # Eq 22.17 TH, problem statement D.1.2

        return self.lambdaET

    def calc_H(self, Rn):
        """
        Calculate Sensible heat flux (D.1.12 in problem statement). The calculated SH flux value will be used
        update Ts (surface temperature) for the current step.

        :param R_n: Net radiation (W m-2).
        :param lambdaET: Evaporation rate (W m-2).

        :return: Sensible heat flux (W m-2).
        """
        self.H = Rn - self.lambdaET  # Given in problem statement D.1.1 & D.1.12

        return self.H

    def calc_Ts(self, Ta_c, ra):
        """
        Calculate surface temperature (Eq. 21.30 in TH). This Ts will be used to calculate LWup in radiation.py script.

        :param Ta_c: Air temperature (deg C).
        :param H: Sensible heat flux (W m-2).
        :param ra: Aerodynamic resistance (s m-1).

        :return: Surface temperature (K).
        """
        Ta_k = Ta_c + 273.17
        self.Ts_k = Ta_k + (self.H * ra / (self.rho_a * self.cp))  # Eq 21.30 TH. Unit in Kelvin

        return self.Ts_k


# This block will help to check the functions if the .py is run as a script (not imported as module)
if __name__ == '__main__':
    # Checking values for grass
    print('=== Grass ===')
    CWB_g = canopy_WB(canopy_type='grass')
    LH = CWB_g.calc_LH(Tc=10.68)
    print(f'LH = {LH} W/m2')

    print(f'delta = {CWB_g.calc_delta(Tc=10.68, e_sat=12.697)} kPa/°C')
    print(f'psy. const = {CWB_g.calc_psy_const(LH)} kPa/°C')
    print(f'Cinterim = {CWB_g.calc_Cinterim(precip=21, Cfinal_prev=0)} mm')

    Cactual = CWB_g.calc_Cactual()
    print(f'Cactual = {Cactual} mm')
    print(f'Dcanopy = {CWB_g.calc_Dcanopy()} mm')

    ra = 14.96
    rs = 244.71
    Rn = 251.65
    VPD = 1.07
    lambda_Ei = CWB_g.calc_lambdaEi(ra=ra, Rn=Rn, VPD_ref=VPD, Cactual=Cactual)
    print(f'Lambda Ei = {lambda_Ei} W/m2')

    Cfinal = CWB_g.calc_Cfinal()
    print(f'Cfinal = {Cfinal} mm')

    lambda_Et = CWB_g.calc_lambdaEt(ra=ra, rs=rs, Rn=Rn, VPD_ref=VPD)
    print(f'Lambda Et = {lambda_Et} W/m2')

    total_ET = CWB_g.calc_lambdaET()
    print(f'Total ET = {total_ET} W/m2')

    print(f'Sensible Heat = {CWB_g.calc_H(Rn=Rn)} W/m2')

    Ts = CWB_g.calc_Ts(Ta_c=10.68, ra=ra)
    print(f'Surface temp. = {Ts} K')


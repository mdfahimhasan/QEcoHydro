# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

from atm import calc_esat, calc_VPD
from ra_rs import calc_ra, calc_rs
from Lnet_Rn import calc_lambda, calc_Rn_RC
from vars_constants import EnvConstants, ForestConstants, GrassConstants


def calc_delta(temp, esat):
    """
    Calculate gradient of the saturation vapor pressure curve from saturation vapor pressure and temperature
    (Eq 2.18 in TH).
    :param esat: saturation vapor pressure  (kPa).
    :param temp: temperature in (deg C).

    :return: slope of the saturation vapor pressure curve (kPa/°C).
    """
    delta = 4098 * esat / (temp + 237.3) ** 2  # Eq 2.18 TH, problem statement D.1.1

    return delta


def calc_LH(temp):
    """
    Calculate latent heat of vaporization from temperature (Eq 2.1 in TH)
    :param temp: temperature (in deg C).

    :return: latent heat of vaporization (lambda) (J/kg).
    """
    LH_vap = (2.501 - 0.002361 * temp) * 10 ** 6  # Eq 2.1 TH, problem statement D.1.1

    return LH_vap


def calc_psy_const(LH_vap):
    """
    Calculate psychrometric constant (Eq 2.25 in TH).
    :param LH_vap: latent heat of vaporization (lambda) (J/kg).

    :return: psychrometric constant (kPa/°C).
    """
    cp = EnvConstants.cp
    P = EnvConstants.P

    psy_const = cp * P / (0.622 * LH_vap)  # Eq. 2.25 TH, problem statement D.1.1

    return psy_const


def calc_ET_PET(temp, e, wind_spd, SW_in, Lnet, FC, WP, Su_0, et_type='et',
                canopy_type='grass', params_ra=None, params_rs=None):
    """
    Calculate ET or PET using Penman-Monteith.

    :param temp: Temperature (deg C).
    :param e: Vapor pressure (kpa).
    :param wind_spd: Wind speed (m/s).
    :param SW_in: Incoming shortwave radiation (W m-2).
    :param Lnet: Net longwave radiation (W/m2).
    :param FC: Field capacity.
    :param WP: Wilting point.
    :param Su_0: Initial storage in the unsaturated zone.
    :param et_type: If 'et' will calculate ET (transpiration from canopy). If 'pet', will calculate PET.
    :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the
                          vars_constants.py file. Default set to 'grass'.
    :param params_ra: A dictionary of params zm, d, z0 to adjust manually for aerodynamic resistance.
                        Default set to None to use these vars from the vars_constants.py.
    :param params_rs: A dictionary of params g0 and gc to adjust manually for stomatal resistance.
                        Default set to None to use these vars from the vars_constants.py.

    :return:
        - E_trans: Transpiration in W/m2.
        - E_trans_mm: Transpiration in mm.
        - ra: Aerodynamic resistance value.
        - rs: Stomatal resistance value.
    """
    # atm params
    e_sat = calc_esat(Tc=temp)  # in kpa
    vpd = calc_VPD(e_sat, e)  # in kpa

    # ra
    if et_type == 'et':
        ra = calc_ra(wind_spd=wind_spd, canopy_type=canopy_type, params=params_ra)  # in s m-1
    else:  # (et_type == 'pet')
        ra = 208 / wind_spd  # in s m-1

    # rs
    if et_type == 'et':
        rs = calc_rs(SW_in=SW_in, vpd=vpd, temp=temp, FC=FC, WP=WP, Su_0=Su_0,
                     canopy_type=canopy_type, params=params_rs)  # in s m-1
    else:  # (et_type == 'pet')
        rs = 70   # in s m-1

    # A (=Rn)
    Rn_RC = calc_Rn_RC(SW_in, Lnet)   # in W m-2

    # calculate other params needed for ET estimation
    delta = calc_delta(temp=temp, esat=e_sat)
    LH_vap = calc_LH(temp=temp)
    gamma = calc_psy_const(LH_vap)
    rho_air = EnvConstants.rho_a
    cp = EnvConstants.cp

    # calculate ET or PET using Penman-Monteith
    E_trans = (delta * Rn_RC + rho_air * cp * vpd / ra) / (delta + gamma * (1 + rs / ra))  # in W m-2

    # conversion factor from J/kg to W/m2
    Wm2_mm = LH_vap ** -1 * 86400

    E_trans_mm = E_trans * Wm2_mm

    return E_trans, E_trans_mm, ra, rs


def SU_eq_mod(Su_0, P, Sumax, alpha, beta, S_canopy_old,
              temp, e, wind_spd, SW_in, Lnet, FC, WP,
              canopy_type='grass', params_ra=None, params_rs=None):
    """
    Unsaturated Zone Reservoir (Su).

    ** This modified version uses Penman-Monteith eqn to calculate ET, instead of a
    a combination of PET and soil moisture availability.

    :param Su_0: Initial storage in the unsaturated zone.
    :param P: Precipitation input.
    :param Sumax: Maximum storage capacity of the unsaturated zone.
    :param alpha: Fraction of precipitation that does not contribute to immediate runoff.
    :param beta: Split between recharge and overland flow.
    :param S_canopy_old: Initial canopy storage (mm).
    :param temp: temperature (deg celsius).
    :param e: vapor pressure.
    :param wind_spd: wind speed.
    :param SW_in: incoming shortwave radiation (W/m2).
    :param Lnet: Net longwave radiation (W/m2).
    :param FC: Field capacity.
    :param WP: Wilting point.
    :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the vars_constants.py.
    :param params_ra: A dictionary of params zm, d, z0 to adjust manually for aerodynamic resistance.
                      Default set to None to use these vars from the vars_constants.py.
    :param params_rs: A dictionary of params g0 and gc to adjust manually for stomatal resistance.
                      Default set to None to use these vars from the vars_constants.py.

    :returns
    - r: Recharge to the saturated zone (contributing to baseflow)
    - SE: Saturation excess overland flow
    - ET_mm: Total actual evapotranspiration
    - E_int_mm: Evaporation from canopy interception
    - E_trans_mm: Transpiration from canopy
    - E_pot_mm : Potential ET
    - S_canopy: Canopy storage
    - S_dt: Updated storage in the unsaturated zone
    - rain_pass: rain not intercepted by canopy
    - ra: aerodynamic resistance
    - rs: stomatal resistance
    """
    if canopy_type == 'grass':
        S_canopy_max = GrassConstants.Scanmax
    elif canopy_type == 'forest':
        S_canopy_max = ForestConstants.Scanmax

    # # Canopy water balance
    # Temporary canopy storage will be old canopy storage plus coming rainfall
    S_canopy_temp = S_canopy_old + P

    # if S_canopy_temp is greater than canopy's maximum capacity, the remaining of S_canopy_temp will be passed as
    # rain_pass (rain not intercepted by canopy). Otherwise, rain_pass will be zero
    if S_canopy_temp > S_canopy_max:
        rain_pass = S_canopy_temp - S_canopy_max
        S_canopy_temp = S_canopy_max
    else:
        rain_pass = 0

    # # Soil water balance
    # Update unsaturated storage with the part of precipitation that was not intercepted by canopy
    Su_0 = Su_0 + (rain_pass * (1 - alpha))

    # Calculate recharge: excess water moves from the unsaturated zone to the saturated zone if storage exceeds capacity
    if Su_0 < Sumax:
        R = 0  # No recharge if storage is below capacity
    else:
        R = Su_0 - Sumax  # Recharge is the excess water above the storage capacity

    # Update unsaturated zone storage by subtracting recharge
    Su_0 = Su_0 - R

    # # Calculate actual evapotranspiration (transpiration) and potential evapotranspiration
    # based on current unsaturated soil moisture storage using Penman-Monteith eqn
    E_trans, E_trans_mm, ra, rs = calc_ET_PET(temp=temp, e=e, wind_spd=wind_spd, SW_in=SW_in, Lnet=Lnet, FC=FC, WP=WP,
                                              Su_0=Su_0, et_type='et', canopy_type=canopy_type,
                                              params_ra=params_ra, params_rs=params_rs)

    E_pot, E_pot_mm, ra_p, rs_p = calc_ET_PET(temp=temp, e=e, wind_spd=wind_spd, SW_in=SW_in, Lnet=Lnet, FC=FC, WP=WP,
                                          Su_0=Su_0, et_type='pet', canopy_type=canopy_type)

    # # updating canopy water balance
    # if PET is greater than canopy temporary storage (canopy old storage + precip or S-canopy_max), evaporation
    # from intercepted water will be equal to canopy temporary storage
    if E_pot_mm > S_canopy_temp:
        excess_ET = E_pot_mm - S_canopy_temp
        E_int_mm = S_canopy_temp

    # if PET is less than canopy temporary storage (canopy old storage + precip or S-canopy_max), evaporation
    # from intercepted water will be equal to potential evaporation
    else:
        excess_ET = 0
        E_int_mm = E_pot_mm

    S_canopy = S_canopy_temp - E_int_mm  # canopy storage
    Su_0 = Su_0 - E_trans_mm  # unsaturated zone storage
    ET_mm = E_int_mm + E_trans_mm  # total ET (from interception and transpiration)

    S_dt = Su_0  # Final storage in the unsaturated zone
    r = R * beta  # Part of the recharge contributing to the baseflow
    SE = R * (1 - beta)  # Part of the recharge contributing to overland flow

    return r, SE, ET_mm, E_int_mm, E_trans_mm, E_pot_mm, S_canopy, S_dt, rain_pass, ra, rs


def SF_eq(S0, P, alpha, Tf, SE):
    """
    Fast Flow Reservoir (Sf).

    :param S0: Initial storage in the fast flow reservoir.
    :param P: Precipitation input.
    :param alpha: Fraction of precipitation that contributes to quick flow.
    :param Tf: Time parameter for quick flow.
    :param SE: Saturation excess overland flow contributing to quick flow.

    :returns
    - Qf: Quick flow (fast runoff)
    - S_dt: Updated storage in the fast flow reservoir
    - IE: Effective infiltration
    """
    # Effective infiltration: part of the precipitation contributing to the fast flow
    IE = (P * alpha)

    # Update storage in the fast flow reservoir with effective infiltration and saturation excess, then subtract outflow
    S_dt = S0 + (P * alpha) + SE - (S0 / Tf)

    # Calculate quick flow based on initial storage and time parameter
    Qf = S0 / Tf

    return Qf, S_dt, IE


def SS_eq(S0, R, Ts):
    """
    Saturated Zone Reservoir (Ss).

    :param S0: Initial storage in the saturated zone.
    :param R: Recharge from the unsaturated zone.
    :param Ts: Time parameter for slow flow.

    :returns
    - Qs: Slow flow (baseflow)
    - S_dt: Updated storage in the saturated zone
    """

    # Update saturated zone storage with recharge, then subtract outflow
    S_dt = S0 + R - (S0 / Ts)

    # Calculate slow flow based on initial storage and time parameter
    Qs = S0 / Ts

    return Qs, S_dt
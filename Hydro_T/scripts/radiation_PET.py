# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

import numpy as np
import matplotlib.pyplot as plt
from Ra_Rso import calc_Ra_Rso


def calc_lambda(Temp):
    """
    Calculate latent heat of vaporization (lambda) at given temperature.

    :param Temp: Temperature (deg celsius).

    :return: Latent heat of vaporization (kJ/kg).
    """
    # Latent Heat of Vaporization (kJ/kg)
    lambda_ = 1000 * (2.501 - 0.002361 * Temp)  # Latent Heat of Vaporization (kJ/kg); TH eq. 2.1

    return lambda_


def calc_Lnet(SW_in, Rso, vp, temp):
    """
    Calculate net longwave radiation.

    :param SW_in: Incoming solar radiation (W/m2).
    :param Rso: Clear sky solar radiation (W/m2).
    :param vp: Vapor pressure (kpa).
    :param temp: Temperature (deg C).

    :returns: net longwave radiation (W/m2).
    """
    # empirical cloud factor
    # SW_in is calculated by TH eq. 5.16 or given.
    # Rso is calculated using dates, lat, elev (look into Ra_Rso.py for detail)
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

    :return: Net solar radiation (W m-2).
    """
    # Net Radiation over a Reference Crop, in W/m2
    # albedo considered 0.23
    albedo = 0.23
    Rn_RC = SW_in * (1 - albedo) + Lnet  # TH eq. 5.28

    return Rn_RC


def calc_PET(u2, Temp, e, Z, SW_in, Lnet):
    """
    Calculates PET.

    :param u2: wind speed (m/s).
    :param Temp: temperature (deg celsius).
    :param e: vapor pressure ().
    :param Z: elevation (m).
    :param SW_in: incoming shortwave radiation (W/m2).
    :param Lnet: Net longwave radiation (W/m2).

    :return: PET (in W/m2 and mm/day) and net radiation (in W/m2 and mm/day).
    """

    # saturated Vapor Pressure (kPa)
    # temp in deg C
    e_sat = 0.6108 * np.exp((17.27 * Temp) / (237.3 + Temp))

    # vapor pressure deficit (kPa)
    vpd = e_sat - e  # Vapor Pressure Deficit in

    # Slope of esat versus temp curve (kPa/C)
    delta = 4098 * e_sat / (237.3 + Temp)**2

    # pressure as function of elevation (kPa)
    press = 101.3 * ((293 - 0.0065 * Z) / 293)**5.26

    # psychrometric Constant (kPa/degC)
    cp = 1.1013  # specific heat at constant pressure for air kJ/(kg.K)
    lambda_ = calc_lambda(Temp)
    gamma = cp * press / (0.622 * lambda_)  # TH eq. 2.25

    # conversion Factor from W/m2 to mm/day (1000 is density of water in kg/m3)
    Wm2_mm = (lambda_ * 1000)**-1 * 86400

    # Rn_RC
    Rn_RC = calc_Rn_RC(SW_in, Lnet)

    # PET eq
    rho_air = 1.23
    ra = 208 / u2  # aerodynamic resistance
    rs = 70  # stomatal/surface resistance
    E_RC = (delta * Rn_RC + rho_air * cp * vpd / ra) / (delta + gamma * (1 + rs / ra))  # TH eq. 22.18

    # mm/day values
    E_RC_mm = E_RC * Wm2_mm
    Rn_RC_mm = Rn_RC * Wm2_mm

    return E_RC, E_RC_mm, Rn_RC, Rn_RC_mm


# This block will help to check the functions if the .py is run as a script (not imported as module)
if __name__ == '__main__':
    # # datasets
    # 1st 3 columns are year, month, day (except Lat_Lon_Area_Z)
    vp_data = '../data/e.csv'  # vapor pressure (kpa)
    Lnet_data = '../data/L_net.csv'  # Net longwave radiation
    precip_data = '../data/Precip.csv'
    Lat_lon_area_Z_data = '../data/Lat_Lon_Area_Z.csv'  # Latitude, longitude, area (m2), and elevation data (m)
    Q_data = '../data/Q_mm.csv'  # streamflow (mm)
    Sin_data = '../data/S_in.csv'  # solar radiation (W/m2)
    Temp_data = '../data/temp.csv'  # temperature (deg C)
    u2_data = '../data/u2.csv'  # wind speed at 2m (m/s)

    # # load datasets
    # removing 1st 3 columns, i.e., year, month, day (except Lat_Lon_Area_Z)
    vp_arr = np.loadtxt(vp_data, delimiter=',')[:, 3:]
    Lnet_arr = np.loadtxt(Lnet_data, delimiter=',')[:, 3:]
    precip_arr = np.loadtxt(precip_data, delimiter=',')[:, 3:]
    Lat_lon_area_Z_arr = np.loadtxt(Lat_lon_area_Z_data, delimiter=',')
    Q_arr = np.loadtxt(Q_data, delimiter=',')[:, 3:]
    Sin_arr = np.loadtxt(Sin_data, delimiter=',')[:, 3:]
    Temp_arr = np.loadtxt(Temp_data, delimiter=',')[:, 3:]
    u2_arr = np.loadtxt(u2_data, delimiter=',')[:, 3:]

    # dates, lat, elev
    dates = np.loadtxt(vp_data, delimiter=',')[:, :3]
    lat = Lat_lon_area_Z_arr[:, 0]
    elev = Lat_lon_area_Z_arr[:, 3]

    # # calculate Ra and Rso
    Ra, Ra_Wm2, Rso, Rso_Wm2 = calc_Ra_Rso(dates=dates, lat=lat, elev=elev)

    # # calculate net longwave
    Lnet = calc_Lnet(SW_in=Sin_arr, Rso=Rso_Wm2, vp=vp_arr, temp=Temp_arr)  # in W/m2

    # # Generate net radiation for reference crop
    # this step is calculated within the calc_PET function in the next step

    # # Generate daily PET
    E_RC, E_RC_mm, Rn_RC, Rn_RC_mm = calc_PET(u2=u2_arr, Temp=Temp_arr, e=vp_arr, Z=elev, SW_in=Sin_arr, Lnet=Lnet)

    # # plots
    # from days 500 to 800
    from_index = 500
    to_index = 800

    # catchment no
    catchment_idx = 12  # for python indexing set to (catchment no - 1)

    # plot PET and net radiation for reference crop
    plt.figure(1)
    plt.clf()
    plt.subplot(1, 2, 1)
    plt.plot(E_RC[from_index:to_index, catchment_idx], '.b', label='Ref. Crop ET')
    plt.plot(Rn_RC[from_index:to_index, catchment_idx], '.-k', label='Net Radiation')
    plt.legend()
    plt.ylabel('W/m2')
    plt.xlabel('days')

    plt.subplot(1, 2, 2)
    plt.plot(E_RC_mm[from_index:to_index, catchment_idx], '.b', label='Ref. Crop ET')
    plt.plot(Rn_RC_mm[from_index:to_index, catchment_idx], '.-k', label='Net Radiation')
    plt.legend()
    plt.ylabel('mm/day')
    plt.xlabel('days')
    plt.suptitle(f'catchment: {catchment_idx + 1}')
    plt.tight_layout()

    plt.savefig(f"../plots/E_RC_vs_RN_RC_catchment_{catchment_idx + 1}.png")
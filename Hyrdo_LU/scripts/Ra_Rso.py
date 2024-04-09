# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

import numpy as np
import matplotlib.pyplot as plt


def calc_Ra_Rso(dates, lat, elev):
    """
    Calculate extraterrestrial radiation (Ra) and clear sky solar radiation (Rso).

    :param dates: numpy array of dates in format [year, month, day].
    :param lat: Latitude in degrees.
    :param elev: Elevation in meters.

    :returns:
            Ra: Extraterrestrial radiation in (MJ m^-2 day^-1).
            Ra_Wm2: Extraterrestrial radiation in (W/m2).
            Rso: Clear sky solar radiation in (MJ m^-2 day^-1).
            Rso_Wm2: Clear sky solar radiation in (W/m2).
    """
    # making Julian dates
    years = np.unique(dates[:, 0])  # unique list of years
    jday_all = np.array([])  # empty array to store julian days

    # there might be leap years in the years.
    # so counting number of days in each year and assigning julian days for the year.
    # finally concatenating to the total list of julian days over all years
    for year in years:
        num_days = len(np.where(dates[:, 0] == year)[0])
        jday_new = np.arange(1, num_days + 1)
        jday_all = np.concatenate((jday_all, jday_new))

    # calculating eccentricity factor (dr), solar declination angle (delta), and Latitude in radians (Lat_rad)
    dr = 1 + 0.033 * np.cos(2 * np.pi * jday_all / 365)  # TH eq. 5.5
    delta = 0.409 * np.sin(2 * np.pi * jday_all / 365 - 1.39)  # TH eq. 5.8
    Lat_rad = np.radians(lat)  # Convert latitude to radians

    # empty array to store Ra or Rso results
    Ra = np.zeros((len(jday_all), len(Lat_rad)))
    Rso = np.zeros((len(jday_all), len(Lat_rad)))

    # calculate Ra or Rso (unit MJ m^-2 day^-1)
    for i, lat in enumerate(Lat_rad):
        ws = np.arccos(-np.tan(delta) * np.tan(lat))  # sunset hour angle; TH eq. 5.12
        Ra[:, i] = (24 * 60 / np.pi) * 0.08202 * dr * (ws * (np.sin(lat) * np.sin(delta)) + np.sin(ws) * (np.cos(lat) * np.cos(delta)))  # TH eq. 5.14
        Rso[:, i] = (0.75 + 2 * 10 ** -5 * elev[i]) * Ra[:, i]

    Rso_Wm2 = Rso / (10 ** -6 * 60 * 60 * 24)  # Convert to W/m2
    Ra_Wm2 = Ra / (10 ** -6 * 60 * 60 * 24)  # Convert to W/m2

    return Ra, Ra_Wm2, Rso, Rso_Wm2


# This block will help to check the functions if the .py is run as a script (not imported as module)
if __name__ == '__main__':
    # datasets
    # 1st 3 columns are year, month, day (except Lat_Lon_Area_Z)
    vp_data = '../data/e.csv'  # vapor pressure (kpa)
    Sin_data = '../data/S_in.csv'  # solar radiation (W/m2)
    Lat_lon_area_Z_data = '../data/Lat_Lon_Area_Z.csv'  # Latitude, longitude, area (m2), and elevation data (m)

    # load datasets
    # removing 1st 3 columns, i.e., year, month, day (except Lat_Lon_Area_Z)
    vp_arr = np.loadtxt(vp_data, delimiter=',')[:, 3:]
    Sin_arr = np.loadtxt(Sin_data, delimiter=',')[:, 3:]
    Lat_lon_area_Z_arr = np.loadtxt(Lat_lon_area_Z_data, delimiter=',')

    # dates, lat, elev
    dates = np.loadtxt(vp_data, delimiter=',')[:, :3]
    lat = Lat_lon_area_Z_arr[:, 0]
    elev = Lat_lon_area_Z_arr[:, 3]

    # calculate Ra, Rso
    Ra, Ra_Wm2, Rso, Rso_Wm2 = calc_Ra_Rso(dates, lat, elev)

    # # plots
    # from days 500 to 800, catchment 8
    from_index = 500
    to_index = 800

    # catchment no
    catchment_idx = 9  # for python indexing set to (catchment no - 1)

    # Plot solar radiations
    plt.figure(1)
    plt.clf()
    plt.plot(Ra_Wm2[from_index:to_index, catchment_idx], '.b', label='Extraterrestrial solar radiation')
    plt.plot(Rso_Wm2[from_index:to_index, catchment_idx], '.k', label='Clear sky solar radiation')
    plt.plot(Sin_arr[from_index:to_index, catchment_idx], '.-r', label='S_in')
    plt.legend()
    plt.ylabel('W/m2')
    plt.xlabel('days')
    plt.suptitle(f'catchment: {catchment_idx + 1}')
    plt.tight_layout()

    plt.savefig(f"../plots/solar_radiation_catchment_{catchment_idx + 1}.png")
import numpy as np
from model import model
from Ra_Rso import calc_Ra_Rso
from reservoirs_ET import calc_ET_PET
from Lnet_Rn import calc_Lnet, calc_Rn_RC
from aggregate import aggregate, make_means_new, make_means

from plots import make_Budyko_plot, plot_P_PET_Temp, plot_Q_monthly, plot_Q_daily, plot_Q_SM_daily

# # datasets
# 1st 3 columns are year, month, day (except Lat_Lon_Area_Z)
# the 15 columns after yy-mm-dd are for 15 different catchments
vp_data = '../data/e.csv'  # vapor pressure (kpa)
precip_data = '../data/Precip.csv'
Lat_lon_area_Z_data = '../data/Lat_Lon_Area_Z.csv'  # Latitude, longitude, area (m2), and elevation data (m)
Q_data = '../data/Q_mm.csv'  # streamflow (mm)
Sin_data = '../data/S_in.csv'  # solar radiation (W/m2)
Temp_data = '../data/temp.csv'  # temperature (deg C)
u2_data = '../data/u2.csv'  # wind speed at 2m (m/s)

# # load datasets
# removing 1st 3 columns, i.e., year, month, day (except Lat_Lon_Area_Z)
vp_arr = np.loadtxt(vp_data, delimiter=',')[:, 3:]
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

# **********************************************************************************************************************
# # calculate Ra and Rso (unit W/m2)
Ra, Ra_Wm2, Rso, Rso_Wm2 = calc_Ra_Rso(dates=dates, lat=lat, elev=elev)
print('Step 1: Estimated extraterrestrial and clear sky radiation radiation...')

# # calculate net longwave
Lnet = calc_Lnet(SW_in=Sin_arr, Rso=Rso_Wm2, vp=vp_arr, temp=Temp_arr, region='arid')  # in W/m2
print('Step 2: Estimated net Longwave radiation...')

# # Generate net radiation for reference crop
# this step is calculated within the calc_PET function in the next step
print('Step 3: Estimated net radiation...')
Rn_RC = calc_Rn_RC(SW_in=Sin_arr, Lnet=Lnet)   # in W/m2


# # calculate PET
canopy_type = 'forest'
print('Step 4: Calculating PET...')
PET, PET_mm, _, _ = calc_ET_PET(temp=Temp_arr, e=vp_arr, wind_spd=u2_arr, SW_in=Sin_arr,
                                Lnet=Lnet, et_type='pet', canopy_type=canopy_type,
                                FC=None, WP=None, Su_0=None)


# # aggregating basin PET, Q, P, and temperature values at monthly and annual scale
PET_mean_monthly_sum, PET_mean_annual = make_means_new(PET_mm, dates)
Q_mean_monthly_sum, Q_mean_annual = make_means_new(Q_arr, dates)
P_mean_monthly_sum, P_mean_annual = make_means_new(precip_arr, dates)
Temp_mean_monthly, Temp_mean_annual = make_means(Temp_arr, dates)

# # Components of Budyko plot
PET_P = PET_mean_annual / P_mean_annual  # aridity index
E_P = (P_mean_annual - Q_mean_annual) / P_mean_annual  # evaporative fraction

# # # Simulation 1
canopy_type = 'grass'
catchment_no = 10
catchment_idx = catchment_no - 1

# # climate inputs
precip_ct = precip_arr[:, catchment_idx]
temp_ct = Temp_arr[:, catchment_idx]
e_ct = vp_arr[:, catchment_idx]
wind_spd_ct = u2_arr[:, catchment_idx]
SW_in_ct = Sin_arr[:, catchment_idx]
Lnet_ct = Lnet[:, catchment_idx]

# # params
mir = 5  # Maximum infiltration rate  [300]
Su_max = 250  # total water capacity [50-300]
Ts = 20  # time parameter for slowflow [20-100]
Tf = 1  # time parameter for quickflow [1-3]
beta = 1  # split between recharge and overland flow.

params_ra = {'zm': 22,     # default
             'd': 14.644,  # default
             'z0': 2}  # Changed
params_rs = {'g0': 5,     # Changed
             'gc': 1}      # default

FC = 0.35 * Su_max
WP = 0.11 * Su_max

# # model run
print('Step 5: Running model simulation...')
output_sim1 = model(P=precip_ct, temp=temp_ct, e=e_ct, wind_spd=wind_spd_ct, SW_in=SW_in_ct, Lnet=Lnet_ct,
                    FC=FC, WP=WP, mir=mir, Su_max=Su_max, Ts=Ts, Tf=Tf, beta=beta, canopy_type=canopy_type,
                    params_ra=params_ra, params_rs=params_rs)

output_sim1 = aggregate(output_dict=output_sim1, dates=dates)

# # Calculating E/P and PET/P for selected catchment
PET_P_SIM1 = output_sim1['PET_mean_annual'] / output_sim1['P_mean_annual']
E_P_SIM1 = (output_sim1['P_mean_annual'] - output_sim1['QT_mean_annual']) / output_sim1['P_mean_annual']

# # Figure 1: Budyko plot
make_Budyko_plot(E_P=E_P, PET_P=PET_P, catchment_no=10,
                 E_P_sim=E_P_SIM1, PET_P_sim=PET_P_SIM1,
                 savepath=f'../plots/BP_catchment_{catchment_no}_sim1.png')

# # Figure 2: P vs PET vs Temp
plot_P_PET_Temp(P_mean_monthly_sum, PET_mean_monthly_sum, Temp_mean_monthly,
                catchment_no=10, save_path=f'../plots/P_PET_Temp_catchment_{catchment_no}_sim1.png')

# # Figure 3: Q plots monthly
plot_Q_monthly(QT_mean_monthly_observed=Q_mean_monthly_sum, QT_mean_monthly_simulated=output_sim1['QT_mean_mon'],
               QS_mean_monthly_simulated=output_sim1['QS_mean_mon'], QF_mean_monthly_simulated=output_sim1['QF_mean_mon'],
               catchment_no=10, save_path=f'../plots/Q_monthly_catchment_{catchment_no}_sim1.png')

# # Figure 4: Q plots daily for a year
plot_Q_daily(QT_daily_observed=Q_arr, QT_daily_simulated=output_sim1['QT'],
             from_day=365, to_day=730, year=1981,
             catchment_no=10, save_path=f'../plots/Q_daily_catchment_{catchment_no}_sim1.png')

# # Figure 5: SM vs Q plots daily in a year

plot_Q_SM_daily(SM_daily_simulated=output_sim1['Su'], QT_daily_simulated=output_sim1['QT'],
                from_day=365, to_day=730, year=1981,
                save_path=f'../plots/Q_vs_SM_daily_catchment_{catchment_no}_sim1.png')

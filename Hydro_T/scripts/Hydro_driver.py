# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

from radiation_PET import *
from toy_model import toymodel

# # datasets
# 1st 3 columns are year, month, day (except Lat_Lon_Area_Z)
# the 15 columns after yy-mm-dd are for 15 different catchments
vp_data = '../data/e.csv'  # vapor pressure (kpa)
precip_data = '../data/Precip.csv'
Lat_lon_area_Z_data = '../data/Lat_Lon_Area_Z.csv'  # Latitude, longitude, area (m2), and elevation data (m)
Q_data = '../data/Q_mm.csv'  # streamflow (mm)
Sin_data = '../data/S_in.csv'  # solar radiation (W/m2)
Temp_data = '../data/Temp.csv'  # temperature (deg C)
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
Lnet = calc_Lnet(Sin=Sin_arr, Rso=Rso_Wm2, vp=vp_arr, temp=Temp_arr)  # in W/m2
print('Step 2: Estimated net Longwave radiation...')

# # Generate net radiation for reference crop
# this step is calculated within the calc_PET function in the next step
print('Step 3: Estimated net radiation...')

# # Generate daily PET (unit available in W/m2 and mm; choose with caution)
E_RC, E_RC_mm, Rn_RC, Rn_RC_mm = calc_PET(u2=u2_arr, Temp=Temp_arr, e=vp_arr, Z=elev, SW_in=Sin_arr, Lnet=Lnet)
print('Step 4: Estimated PET...')

# # Run hydrologic model for ""selected catchment""

# catchment no
catchment_idx = 12  # for python indexing set to (catchment no - 1)

# preparing data for catchment
precip_catchment = precip_arr[:, catchment_idx]
E_RC_catchment = E_RC_mm[:, catchment_idx]

# the reservoir operations will be implemented inside the hydrologic model
Ea, QF, R, QS, QT, Sf, Su, Ss, St, AL, IE, SE = \
    toymodel(P=precip_catchment, Ep=E_RC_catchment, mir=50, Su_max=450, Ts=5, Tf=1, beta=1)
print(f'Step 5: Ea, QF, R, QS, QT, Sf, Su, Ss, St, AL, IE, SE for catchment {catchment_idx + 1}...')

# **********************************************************************************************************************

# # plot
# from days 500 to 800
from_index = 500
to_index = 800

fig, ax = plt.subplots(1, 2, figsize=(25, 8))
plt.rcParams.update({'font.size': 20})

ax[0].plot(Q_arr[from_index:to_index, catchment_idx], '-b', label='Observed')
ax[0].plot(QT[from_index:to_index], '-k', label='Simulated')
ax[0].set_ylabel('mm', fontsize=20)
ax[0].set_xlabel('days', fontsize=20)
ax[0].tick_params(axis='both', labelsize=20)
ax[0].legend()

Su_max = 450  # maximum soil moisture holding capacity
h1, = ax[1].plot(Su[from_index:to_index] / Su_max, '.-b', label='Su/Su_max')
ax[1].set_xlabel('days', fontsize=20)
ax[1].set_ylabel('Su/Su_max \n [-]', color='b', fontsize=20)  # Set y-label color to blue to match the data
ax[1].tick_params(axis='y', colors='b')  # Set tick color to blue to match the data
ax[1].tick_params(axis='both', labelsize=20)
ax[1] = plt.twinx()
h2, =ax[1].plot(Ea[from_index:to_index], '.-r', label='ET')
ax[1].set_ylabel('mm', color='r', fontsize=20)  # Set y-label color to red to match the data
ax[1].tick_params(axis='y', colors='r')  # Set tick color to red to match the data
ax[1].tick_params(axis='both', labelsize=20)

handles = [h1, h2]
labels = [h.get_label() for h in handles]
ax[1].legend(handles, labels)

fig.suptitle(f'catchment: {catchment_idx + 1}')

plt.savefig(f"../plots/Q_ET_catchment_{catchment_idx + 1}.png")














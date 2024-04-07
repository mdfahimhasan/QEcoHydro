# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

import os
import pandas as pd

from atm import atmosphere
from aeroDyn_res import aerodynamic_res
from canopy_res import surface_res
from radiation import radiation
from canopy_WB import canopy_WB
from soil_WB import soil_WB


def run_SVAT(forcing_csv, canopy_type, output_csv):
    """
    Run SVAT model with the given forcing csv for particular canopy type.

    :param forcing_csv: Filepath of forcing data csv.
    :param canopy_type: Canopy type from 'grass' or 'forest'. Currently, the model can't handle any other canopy type.
    :param output_csv: Output csv filepath.

    :return: A csv with model results.
    """
    # reading forcing data
    data_df = pd.read_csv(forcing_csv)

    # # # Part A: Aerodynamic Parameterization
    ar = aerodynamic_res(canopy_type=canopy_type)  # initializing aerodynamic_res class instance
    d = ar.calc_d()  # zero plane displacement height (m)
    z0 = ar.calc_z0()  # aerodynamic roughness of crop (m)
    data_df['ra'] = data_df['wind_speed'].apply(ar.calc_ra)

    # # # Part B: Surface Resistance Parameterization
    atm = atmosphere()  # initializing atm class instance
    data_df['e'] = data_df['q'].apply(atm.calc_e)  # vapor pressure (kPa)
    data_df['e_sat'] = data_df['Ta'].apply(atm.calc_esat)  # saturation vapor pressure (kPa)
    data_df['VPD'] = data_df.apply(lambda row: atm.calc_VPD(row['e_sat'], row['e']), axis=1)  # vapor pressure deficit (kPa)

    sr = surface_res(canopy_type=canopy_type)  # initializing surface_res class instance
    data_df['gR'] = data_df['SR_down'].apply(sr.calc_gR)  # resistance to soil heat flux (none)
    data_df['gD'] = data_df['VPD'].apply(sr.calc_gD)  # resistance to vapor pressure deficit (none)
    data_df['gT'] = data_df['Ta'].apply(sr.calc_gT)  # resistance to temperature (none)

    # # gSM has to estimated at each step as SM is dynamically estimated at every step using canopy & soil moisture water balance
    # # gS has to be calculated at each step as it depends of dynamic gSM

    # # # Part C: Radiation Parameterization
    rad = radiation(canopy_type=canopy_type)  # initializing radiation class instance

    # # both LW_up and Rn has to calculated at each step as they depend on Ts (surface temperature)
    # # Ts is dependent of sensible heat (H) at each step

    # # # Part D: Canopy Water Balance Parameterization (Includes Interception)
    cwb = canopy_WB(canopy_type=canopy_type)  # initializing canopy_WB class instance
    data_df['Lambda'] = data_df['Ta'].apply(cwb.calc_LH)  # latent heat of vaporization (lambda) (J/kg)
    data_df['delta'] = data_df.apply(lambda row: cwb.calc_delta(row['Ta'], row['e_sat']),
                                     axis=1)  # slope of the saturation vapor pressure curve (kPa/°C)
    data_df['psy_const'] = data_df['Lambda'].apply(cwb.calc_psy_const)  # psychrometric constant (kPa/°C)

    # # Cinterim, Cactual, Dcanopy, lambdaEi, Cfinal, lambdaEt, lambdaET, H, and Ts have to estimated at each step dynamically

    # # # Part E: Soil Water Balance Parameterization
    swb = soil_WB(canopy_type=canopy_type)  # initializing soil_WB class instance

    # # SMnew have to be estimated at each time step as it depends on fluxes estimated at the canopy_WB class

    # # # # running loop for each time step to calculate step-dependent dynamic values
    for t_step in range(len(data_df)):
        print(f'running model for {t_step=}, hr = {data_df.at[t_step, "Total Hour"]}')

        # # # Part B
        # # gSM, gS, and rs
        if t_step == 0:
            data_df.at[t_step, 'gSM'] = sr.calc_gSM(SM=sr.SMinit)  # SM set to SMinit at the 1st time step
        else:
            data_df.at[t_step, 'gSM'] = sr.calc_gSM(SM=data_df.at[(t_step-1), 'SMnew'])   # from 2nd step SM set to SM_last from previous time step

        # extracting gD, gT, and gSM values for this time step and setting to the sr class for method execution
        sr.gR = data_df.at[t_step, 'gR']
        sr.gD = data_df.at[t_step, 'gD']
        sr.gT = data_df.at[t_step, 'gT']
        sr.gSM = data_df.at[t_step, 'gSM']

        data_df.at[t_step, 'gS'] = sr.calc_gS()  # total surface conductance (mm s-1)
        data_df.at[t_step, 'rs'] = sr.calc_rs(gS=data_df.at[t_step, 'gS'])  # total surface resistance (s m-1)

        # # # Part C
        # # Lu and Rn
        if (t_step == 0) | (t_step == 1):
            data_df.at[t_step, 'LW_up'] = rad.LW_up_init  # for 1st two steps LW_up_init set to initial values
        else:
            Ts1_k = data_df.at[(t_step - 1), 'Ts']  # surface temp of 1 step before (K)
            Ts2_k = data_df.at[(t_step - 2), 'Ts']  # surface temp of 2 step before (K)
            data_df.at[t_step, 'LW_up'] = rad.calc_LW_up(Ts1_k, Ts2_k)  # longwave radiation from surface (W m-2)

        data_df.at[t_step, 'Rn'] = rad.calc_Rn(SW_in=data_df.at[t_step, 'SR_down'],
                                               LW_up=data_df.at[t_step, 'LW_up'],
                                               LW_down=data_df.at[t_step, 'LW_down'])  # net solar radiation (W m-2)

        # # # Part D
        # # Cinterim, Cactual, Dcanopy, lambdaEi, Cfinal, lambdaEt, lambdaEt, H, Ts
        if t_step == 0:
            data_df.at[t_step, 'Cinterim'] = 0
        else:
            data_df.at[t_step, 'Cinterim'] = cwb.calc_Cinterim(precip=data_df.at[t_step, 'precip'],
                                                               Cfinal_prev=data_df.at[(t_step-1), 'Cfinal'])  # interim canopy storage (mm)

        # extracting Cinterim, delta, Lambda, psy_const values for this time step and
        # setting them in the cwb class for method execution
        cwb.Cinterim = data_df.at[t_step, 'Cinterim']
        cwb.delta = data_df.at[t_step, 'delta']
        cwb.LH = data_df.at[t_step, 'Lambda']
        cwb.psy_const = data_df.at[t_step, 'psy_const']

        data_df.at[t_step, 'Cactual'] = cwb.calc_Cactual()  # actual canopy storage (mm)
        data_df.at[t_step, 'Dcanopy'] = cwb.calc_Dcanopy()  # canopy drainage (mm)
        data_df.at[t_step, 'lambdaEi'] = cwb.calc_lambdaEi(ra=data_df.at[t_step, 'ra'],
                                                           Rn=data_df.at[t_step, 'Rn'],
                                                           VPD_ref=data_df.at[t_step, 'VPD'],
                                                           Cactual=data_df.at[t_step, 'Cactual'])  # evaporation of intercepted rainfall with free water canopy (W m-2)
        data_df.at[t_step, 'Cfinal'] = cwb.calc_Cfinal()  # final canopy storage (mm)
        data_df.at[t_step, 'lambdaEt'] = cwb.calc_lambdaEt(ra=data_df.at[t_step, 'ra'],
                                                           rs=data_df.at[t_step, 'rs'],
                                                           Rn=data_df.at[t_step, 'Rn'],
                                                           VPD_ref=data_df.at[t_step, 'VPD'])  # transpiration from canopy (W m-2)
        data_df.at[t_step, 'lambdaET'] = cwb.calc_lambdaET()  # total ET (W m-2)
        data_df.at[t_step, 'H'] = cwb.calc_H(Rn=data_df.at[t_step, 'Rn'])  # sensible heat flux (W m-2)
        data_df.at[t_step, 'Ts'] = cwb.calc_Ts(Ta_c=data_df.at[t_step, 'Ta'],
                                               ra=data_df.at[t_step, 'ra'])  # surface temperature (K)

        # # # Part E
        # #  SMnew and SMlast
        if t_step == 0:
            data_df.at[t_step, 'SMlast'] = sr.SMinit
        else:
            data_df.at[t_step, 'SMlast'] = data_df.at[(t_step-1), 'SMnew']

        data_df.at[t_step, 'SMnew'] = swb.calc_SMnew(Dcanopy=data_df.at[t_step, 'Dcanopy'],
                                                     Cactual=data_df.at[t_step, 'Cactual'],
                                                     lambdaET=data_df.at[t_step, 'lambdaET'],
                                                     LH=cwb.LH,
                                                     SMlast=data_df.at[t_step, 'SMlast'])

    # saving dataframe as csv
    data_df.to_csv(output_csv, index=False)


if __name__ == '__main__':
    forcing_csv = '../forcing_data/forcings.csv'
    canopy_type = 'forest'
    output_csv = '../results/forest.csv'
    outdir = os.path.dirname(output_csv)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    run_SVAT(forcing_csv, canopy_type, output_csv)



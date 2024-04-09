import numpy as np
from reservoirs_ET import SF_eq, SU_eq_mod, SS_eq


def model(P, temp, e, wind_spd, SW_in, Lnet, FC, WP, mir, Su_max, Ts, Tf, beta,
          canopy_type='grass', params_ra=None, params_rs=None):
    """
    Model for Hydrologic Simulation.

    :param P: (list) Precipitation (mm).
    :param temp: temperature (deg celsius).
    :param e: vapor pressure.
    :param wind_spd: wind speed.
    :param SW_in: incoming shortwave radiation (W/m2).
    :param Lnet: Net longwave radiation (W/m2).
    :param FC: Field capacity.
    :param WP: Wilting point.
    :param mir: Maximum infiltration rate.
    :param Su_max: Max storage capacity of unsaturated zone.
    :param Ts: Time parameter for slow flow.
    :param Tf: Time parameter for quick flow.
    :param beta: Split between recharge and overland flow.
    :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the vars_constants.py file.
    :param params_ra: A dictionary of params zm, d, z0 to adjust manually for aerodynamic resistance.
                      Default set to None to use these vars from the vars_constants.py.
    :param params_rs: A dictionary of params g0 and gc to adjust manually for stomatal resistance.
                      Default set to None to use these vars from the vars_constants.py.

    :returns
    - Ea: Actual evapotranspiration
    - QF: Quick flow (fast runoff)
    - R: Recharge to the saturated zone (contributing to baseflow)
    - QS: Slow flow (baseflow)
    - QT: Total flow (QF + QS)
    - Sf: Storage in the fast flow reservoir
    - Su: Storage in the unsaturated zone
    - Ss: Storage in the saturated zone
    - St: Total storage (Su + Ss)
    - alpha: Fraction of precipitation that contributes to runoff
    - IE: Effective infiltration
    - SE: Saturation excess overland flow
    """
    M = len(P)  # number of data points/iterations/steps

    # Initialize output arrays
    QF, QS, QT, R, Sf, Su, Ss, Ea, alpha, IE, SE, E_int, E_trans, E_pot, S_canopy, pot_inf = \
        [np.zeros(M) for _ in range(16)]  # the integer depends on the number of variable stores we want to open

    ET_vars = {'ra': np.zeros(M), 'rs': np.zeros(M)}

    # Initial conditions
    S0 = 200  # Initial storage for all reservoirs

    # Initial conditions for storage in each reservoir
    Su_dt = S0
    Ss_dt = S0
    Sf_dt = S0
    S_canopy_old = 0

    # Main flow routine
    for t in range(M):

        # Fraction of precipitation contributing to runoff
        if P[t] > 0:
            alpha[t] = 1 - (1 - np.exp(-P[t] / mir)) / (P[t] / mir)
        else:
            alpha[t] = 0

        # Update components of unsaturated zone reservoir (recharge, evapotranspiration, unsaturated storage, and saturation excess)
        r, se, ET_mm, E_int_mm, E_trans_mm, E_pot_mm, Scanopy, S_dt, rain_pass, ra, rs = \
            SU_eq_mod(Su_0=Su_dt, P=P[t], Sumax=Su_max, alpha=alpha[t], beta=beta, S_canopy_old=S_canopy_old,
                      temp=temp[t], e=e[t], wind_spd=wind_spd[t], SW_in=SW_in[t], Lnet=Lnet[t], FC=FC, WP=WP,
                      canopy_type=canopy_type, params_ra=params_ra, params_rs=params_rs)

        R[t], Ea[t], Su[t], SE[t], E_int[t], E_trans[t], E_pot[t], S_canopy[t], S_canopy_old, pot_inf[t] = \
            r, ET_mm, S_dt, se, E_int_mm, E_trans_mm, E_pot_mm, Scanopy, Scanopy, rain_pass

        ET_vars['ra'][t] = ra
        ET_vars['rs'][t] = rs

        # Update saturated zone reservoir (slow flow and saturated zone storage)
        Qs, Ss_dt = SS_eq(Ss_dt, r, Ts)

        QS[t], Ss[t] = Qs, Ss_dt

        # Update components of fast flow reservoir (quick flow, fast flow reservoir storage, and effective infiltration)
        Qf, Sf_dt, ie = SF_eq(Sf_dt, P[t], alpha[t], Tf, SE[t])

        QF[t], Sf[t], IE[t] = Qf, Sf_dt, ie

    QT = QF + QS  # Total flow is the sum of quick flow and slow flow
    St = Su + Ss  # Total storage is the sum of unsaturated and saturated zone storages

    output_dict = {
                    'Ea': Ea, 'QF': QF, 'R': R, 'QS': QS, 'QT': QT,
                    'Sf': Sf, 'Su': Su, 'Ss': Ss, 'St': St, 'AL': alpha, 'IE': IE,
                    'SE': SE, 'Ei': E_int, 'Et': E_trans, 'Ep': E_pot, 'S_canopy': S_canopy,
                    'pot_inf': pot_inf, 'ET_vars': ET_vars, 'P': P
                  }

    return output_dict
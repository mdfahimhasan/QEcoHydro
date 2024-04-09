# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

import numpy as np
from reservoirs import SF_eq, SU_eq, SS_eq


def toymodel(P, Ep, mir, Su_max, Ts, Tf, beta):
    """
    Toy Model Hydrologic Simulation.

    :param P: (list) Precipitation (mm).
    :param Ep: (list) Potential evapotranspiration (mm).
    :param mir: Maximum infiltration rate.
    :param Su_max: Max storage capacity of unsaturated zone.
    :param Ts: Time parameter for slow flow.
    :param Tf: Time parameter for quick flow.
    :param beta: Split between recharge and overland flow.

    :returns
    - Ea: Actual evapotranspiration
    - QF: Quick flow (fast runoff)
    - r: Recharge to the saturated zone (contributing to baseflow)
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
    QF, QS, SE, QT, r, Sf, Su, Ss, Ea, alpha, IE = [np.zeros(M) for _ in range(11)]

    # Initial conditions
    S0 = 0.2  # Initial storage for all reservoirs

    # Initial conditions for storage in each reservoir
    Su_dt = S0
    Ss_dt = S0
    Sf_dt = S0  # Assume an initial storage value

    # Main flow routine
    for t in range(M):

        # Fraction of precipitation contributing to runoff
        if P[t] > 0:
            alpha[t] = 1 - (1 - np.exp(-P[t] / mir)) / (P[t] / mir)
        else:
            alpha[t] = 0

        # Update components of unsaturated zone reservoir (recharge, evapotranspiration, unsaturated storage, and saturation excess)
        r, se, E, Su_dt = SU_eq(Su_0=Su_dt, P=P[t], Ep=Ep[t], Sumax=Su_max, alpha=alpha[t], beta=beta)

        r[t], Ea[t], Su[t], SE[t] = r, E, Su_dt, se

        # Update saturated zone reservoir (slow flow and saturated zone storage)
        Qs, Ss_dt = SS_eq(Ss_dt, r, Ts)

        QS[t], Ss[t] = Qs, Ss_dt

        # Update components of fast flow reservoir (quick flow, fast flow reservoir storage, and effective infiltration)
        Qf, Sf_dt, ie = SF_eq(Sf_dt, P[t], alpha[t], Tf, SE[t])

        QF[t], Sf[t], IE[t] = Qf, Sf_dt, ie

    QT = QF + QS  # Total flow is the sum of quick flow and slow flow
    St = Su + Ss  # Total storage is the sum of unsaturated and saturated zone storages

    return Ea, QF, r, QS, QT, Sf, Su, Ss, St, alpha, IE, SE

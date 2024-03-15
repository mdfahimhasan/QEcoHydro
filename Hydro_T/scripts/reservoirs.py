import numpy as np


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
    - ie: Effective infiltration
    """
    # Effective infiltration: part of the precipitation contributing to the fast flow
    ie = (P * alpha)

    # Update storage in the fast flow reservoir with effective infiltration and saturation excess, then subtract outflow
    S_dt = S0 + (P * alpha) + SE - (S0 / Tf)

    # Calculate quick flow based on initial storage and time parameter
    Qf = S0 / Tf

    return Qf, S_dt, ie


def SU_eq(Su_0, P, Ep, Sumax, alpha, beta):
    """
    Unsaturated Zone Reservoir (Su).

    :param Su_0: Initial storage in the unsaturated zone.
    :param P: Precipitation input.
    :param Ep: Potential evapotranspiration.
    :param Sumax: Maximum storage capacity of the unsaturated zone.
    :param alpha: Fraction of precipitation that does not contribute to immediate runoff.
    :param beta: Split between recharge and overland flow.

    :returns
    - r: Recharge to the saturated zone
    - se: Saturation excess overland flow
    - E: Actual evapotranspiration
    - S_dt: Updated storage in the unsaturated zone
    """
    # Update unsaturated storage with the part of precipitation not contributing to immediate runoff
    Su_0 = Su_0 + (P * (1 - alpha))

    # Calculate recharge: excess water moves from the unsaturated zone to the saturated zone if storage exceeds capacity
    if Su_0 < Sumax:
        R = 0  # No recharge if storage is below capacity
    else:
        R = Su_0 - Sumax  # Recharge is the excess water above the storage capacity

    # Update unsaturated zone storage by subtracting recharge
    Su_0 = Su_0 - R

    # Calculate actual evapotranspiration based on current storage and potential evapotranspiration
    E = Ep * (Su_0 / Sumax)

    if E > Ep:
        E = Ep  # Ensure actual evapotranspiration does not exceed potential

    # Update unsaturated zone storage by subtracting actual evapotranspiration
    Su_0 = Su_0 - E

    S_dt = Su_0  # Final storage in the unsaturated zone
    r = R * beta  # Part of the recharge contributing to the baseflow
    se = R * (1 - beta)  # Part of the recharge contributing to overland flow
    return r, se, E, S_dt


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
# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

# ===================== Independent Constants =====================

class EnvConstants:
    """Environmental constants for atmospheric and crop conditions with units and descriptions.

    Attributes:
    :param P (float): Atmospheric pressure (kPa).
    :param k (float): von Karman's constant (none).
    :param sigma (float): Stefan-Boltzmann constant (W m-2 K-4).
    :param KR (int): Plant-specific parameter in Eq (24.2) (W m-2).
    :param KD1 (float): Plant-specific parameter in Eq (24.6) (kPa-1).
    :param KD2 (float): Plant-specific parameter in Eq (24.6) (kPa-2).
    :param TL (float): Temperature parameter in Eq (24.4) and Eq (24.5) (K).
    :param T0 (float): Temperature parameter in Eq (24.4) and Eq (24.5) (K).
    :param TH (float): Temperature parameter in Eq (24.4) and Eq (24.5) (K).
    :param KM2 (float): Parameter in Eq (24.6) (mm-1).
    :param rho_a (float): Moist Air density (kg m-3).
    :param cp (float): Specific heat capacity of air at constant pressure (J kg-1 K-1).
    :param gc (float): Canopy cover factor (none).
    :param aT (float): Parameter in temperature stress factor (none).
    """
    P: float = 101.20  # Atmospheric pressure
    k: float = 0.4  # von Karman's constant
    sigma: float = 5.67E-08  # Stefan-Boltzmann constant (W m-2 K-4)
    KR: int = 200  # Parameter in Eq (24.2)
    KD1: float = -0.307  # Parameter in Eq (24.6)
    KD2: float = 0.019  # Parameter in Eq (24.6)
    TL: float = 273.00  # Parameter in Eq (24.4) and Eq (24.5)
    T0: float = 293.00  # Parameter in Eq (24.4) and Eq (24.5)
    TH: float = 313.00  # Parameter in Eq (24.4) and Eq (24.5)
    KM2: float = -0.10  # Parameter in Eq (24.6)
    rho_a: float = 1.23  # Moist Air density
    cp: float = 1013.00  # Specific heat capacity of air at constant pressure
    gc: float = 1.00  # Canopy cover factor
    aT: float = 1.00  # Parameter in temperature stress factor


# ===================== Forest Dependent Constants =====================

class ForestConstants:
    """Forest environmental constants with units and descriptions.

    Attributes:
    :param LAI (float): Leaf Area Index (none).
    :param h (float): Canopy Height (m).
    :param a (float): Albedo (none).
    :param g0 (float): Canopy Specific Constant (mm s-1).
    :param E_surface (float): Emissivity of forest (none).
    :param KM1 (float): Parameter in Eq (24.6) (none).
    :param SMo (float): Maximum soil moisture accessible to roots. Parameter in Eq (24.6) (mm).
    :param SMinit (float): Initial Root-accessible soil moisture (mm).
    :param S (float): Maximum Canopy Water Storage (mm).
    :param zm (float): Measurement height for wind speed (m). Considered equal to measurement height for temperature.
    """
    LAI: float = 4.00  # Leaf Area Index
    h: float = 20.00  # Canopy Height
    a: float = 0.12  # Albedo
    g0: float = 15.00  # Canopy Specific Constant
    E_surface: float = 0.95  # Emissivity of crop
    KM1: float = 3.36E-04  # Parameter in Eq (24.6)
    SMo: float = 80.00  # Maximum soil moisture accessible to roots
    SMinit: float = 40.00  # Initial Root-accessible soil moisture
    S: float = 4.00  # Maximum Canopy Water Storage
    zm: float = 22  # Measurement height for wind speed/temperature (m)
    LW_up_init: float = -350.00  # Initial upward longwave radiation (W m-2)


# ===================== Grass Dependent Constants =====================

class GrassConstants:
    """Grass environmental constants with units and descriptions.

    Attributes:
    :param LAI (float): Leaf Area Index (none).
    :param h (float): Canopy Height (m).
    :param a (float): Albedo (none).
    :param g0 (float): Canopy Specific Constant (mm s-1).
    :param E_surface (float): Emissivity of grass (none).
    :param KM1 (float): Parameter in Eq (24.6) (none).
    :param SMo (float): Maximum soil moisture accessible to roots. Parameter in Eq (24.6) (mm).
    :param SMinit (float): Initial Root-accessible soil moisture (mm).
    :param S (float): Maximum Canopy Water Storage (mm).
    :param zm (float): Measurement height for wind speed (m). Considered equal to measurement height for temperature.
    """
    LAI: float = 2.00  # Leaf Area Index
    h: float = 0.12  # Canopy Height
    a: float = 0.23  # Albedo
    g0: float = 30.00  # Canopy Specific Constant
    E_surface: float = 0.95  # Emissivity of crop
    KM1: float = 0.01865  # Parameter in Eq (24.6)
    SMo: float = 40.00  # Maximum soil moisture accessible to roots
    SMinit: float = 10.00  # Initial Root-accessible soil moisture
    S: float = 2.00  # Maximum Canopy Water Storage
    zm: float = 2  # Measurement height for wind speed/temperature (m)
    LW_up_init: float = -319.00  # Initial upward longwave radiation (W m-2)

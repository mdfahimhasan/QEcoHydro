# Author: Md Fahim Hasan
# Email: Fahim.Hasan@colostate.edu

from vars_constants import ForestConstants, GrassConstants


# Part E: Soil Water Balance Parameterization

class soil_WB:

    def __init__(self, canopy_type='grass'):
        """
        Calculates new value of soil moisture.

        :param SMo: Maximum soil moisture accessible to roots (mm). Crop constant.
        :param S: Maximum Canopy Water Storage (mm).
        :param canopy_type: 'forest' or 'grass'. Will use the text to get values of h, LAI, and zm from the
                                 vars.py file. Default set to 'grass'.

        :param SMnew: New soil moisture (mm). Default set to None.
        """
        # canopy constants
        if canopy_type == 'grass':
            self.SM0 = GrassConstants.SMo
            self.S = GrassConstants.S

        elif canopy_type == 'forest':
            self.SM0 = ForestConstants.SMo
            self.S = ForestConstants.S

        # variables estimated using methods (functions) within the class
        self.SMnew = None

    def calc_SMnew(self, Dcanopy, Cactual, lambdaET, LH, SMlast):
        """
        Calculate new soil moisture (E.1 & E.2 in problem statement)

        :param Dcanopy: Canopy water storage change (mm).
        :param Cactual: Actual canopy water storage (mm).
        :param LH: Latent heat of vaporization (lambda) (J/kg).
        :param lambdaET: Evapotranspiration (mm).
        :param SMlast: Previous soil moisture (mm).

        :return: New soil moisture (mm).
        """
        m_to_mm = 1000
        rho_w = 1000  # Density of water in kg/m^3
        lambdaET_mm = lambdaET * 3600 / LH / rho_w * m_to_mm  # Convert from W m-2 to mm

        # the 2nd part of the equation with the (1-C/s) term is the portion coming through (from soil moisture) that
        # is used for canopy transpiration. As this part is lost to transpiration, its subtracted here but is used a
        # positive term in ET equation (calc_lambdaET() func in canopy_WB.py)
        SMnew = Dcanopy - (1 - Cactual / self.S) * lambdaET_mm  # problem statement E.1

        # Ensure SMnew does not exceed soil moisture holding capacity
        SMnew = SMnew + SMlast

        SMnew = min(SMnew, self.SM0)  # problem statement E.2, make sure its not greater than moisture holding capacity
        self.SMnew = max(SMnew, 0)  # Ensure SMnew is not negative

        return self.SMnew

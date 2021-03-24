import numpy as np

class Constants:
    """ An object defining constants
    """ 

    def __init__(self):

        # Gravitational Constant
        self.G = 4*np.pi**2 #AU^3 yr^-2 M_sun^-1

        # Astronomical Unit in meters
        self.AU = 1.496e+11 # m

        # Seconds in a year
        self.yr = 31557600 #s

        # Solar mass in kg
        self.M_sun = 1.989e30 #kg

        # Earth mass in kg
        self.M_earth = 3.0404327497692654e-06*self.M_sun #kg

        # Days per year
        self.days_per_year = 365.2568983263281 #days


import numpy as np
from bin.classes.AstroLi import CelestialBody, KeplerOrbit

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

class Planets:
    """ An object containing planet constants

    Example:
    p = Planets()
    p.Earth.mass
    >> returns earths mass

    """

    def __init__(self):

        MERCURY_ORBIT = KeplerOrbit(0.38709893, 0.20563069, 7.00487, O=48.33167, w_tilde=77.45645)
        VENUS_ORBIT = KeplerOrbit(0.72333199, 0.00677323, 3.39471, O=76.68069, w_tilde=131.53298)
        EARTH_ORBIT = KeplerOrbit(1.00000011, 0.01671022, 0.00005, O=-11.26064, w_tilde=102.94719)
        MARS_ORBIT = KeplerOrbit(1.52366231, 0.09341233, 1.85061, O=49.57854, w_tilde=336.04084)
        JUPITER_ORBIT = KeplerOrbit(5.20336301, 0.04839266, 1.30530, O=100.55615, w_tilde=14.75385)
        SATURN_ORBIT = KeplerOrbit(9.53707032, 0.05415060, 2.48446, O=113.71504, w_tilde=92.43194)
        URANUS_ORBIT = KeplerOrbit(19.19126393, 0.04716771, 0.76986, O=74.22988, w_tilde=170.96424)
        NEPTUNE_ORBIT = KeplerOrbit(30.06896348, 0.00858587, 1.76917, O=131.72169, w_tilde=44.97135)

        # Masses of planets given in 1/M_Sun
        MERCURY_MASS = 6023600
        VENUS_MASS = 408523.71
        EARTH_MASS = 332946.050895
        MARS_MASS =  3098708
        JUPITER_MASS = 1047.3486
        SATURN_MASS =  3497.898
        URANUS_MASS =  22902.98
        NEPTUNE_MASS =  19412.24
        PLUTO_MASS = 1.352e8
        MOON_MASS = 27068700.387534
        SUN_MASS = 1.0

        self.MERCURY =   CelestialBody(name="Mercury",   typ="Planet", mass=1/MERCURY_MASS,   k_orbit=MERCURY_ORBIT)
        self.VENUS =     CelestialBody(name="Venus",     typ="Planet", mass=1/VENUS_MASS,   k_orbit=VENUS_ORBIT)
        self.EARTH =     CelestialBody(name="Earth",     typ="Planet", mass=(1/EARTH_MASS + 1/MOON_MASS),   k_orbit=EARTH_ORBIT)
        self.MARS =      CelestialBody(name="Mars",      typ="Planet", mass=1/MARS_MASS,   k_orbit=MARS_ORBIT)
        self.JUPITER =   CelestialBody(name="Jupiter",   typ="Planet", mass=1/JUPITER_MASS,   k_orbit=JUPITER_ORBIT)
        self.SATURN =    CelestialBody(name="Saturn",    typ="Planet", mass=1/SATURN_MASS,   k_orbit=SATURN_ORBIT)
        self.URANUS =    CelestialBody(name="Uranus",    typ="Planet", mass=1/URANUS_MASS,   k_orbit=URANUS_ORBIT)
        self.NEPTUNE =   CelestialBody(name="Neptune",   typ="Planet", mass=1/NEPTUNE_MASS,   k_orbit=NEPTUNE_ORBIT)
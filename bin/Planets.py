import numpy as np
from bin.classes.AstroLi import CelestialBody, KeplerOrbit


class Planets:
    """ An object containing planet constants

    Example:
    p = Planets()
    p.Earth.mass
    >> returns earths mass

    """

    def __init__(self):

        MERCURY_ORBIT = KeplerOrbit(0.38709893, 0.20563069, 7.00487, O=48.33167, w_tilde=77.45645, L=252.25084)
        VENUS_ORBIT = KeplerOrbit(0.72333199, 0.00677323, 3.39471, O=76.68069, w_tilde=131.53298, L=181.97973)
        EARTH_ORBIT = KeplerOrbit(1.00000011, 0.01671022, 0.00005, O=-11.26064, w_tilde=102.94719, L=100.46435)
        MARS_ORBIT = KeplerOrbit(1.52366231, 0.09341233, 1.85061, O=49.57854, w_tilde=336.04084, L=355.45332)
        JUPITER_ORBIT = KeplerOrbit(5.20336301, 0.04839266, 1.30530, O=100.55615, w_tilde=14.75385, L=34.40438)
        SATURN_ORBIT = KeplerOrbit(9.53707032, 0.05415060, 2.48446, O=113.71504, w_tilde=92.43194, L=49.94432)
        URANUS_ORBIT = KeplerOrbit(19.19126393, 0.04716771, 0.76986, O=74.22988, w_tilde=170.96424, L=313.23218)
        NEPTUNE_ORBIT = KeplerOrbit(30.06896348, 0.00858587, 1.76917, O=131.72169, w_tilde=44.97135, L=304.88003)
        PLUTO_ORBIT = KeplerOrbit(39.48168677, 0.24880766, 17.14175, O=110.30347, w_tilde=224.06676, L=238.92881)
        SUN_ORBIT = KeplerOrbit(0, 0, 0, 0, 0, 0)

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
        self.PLUTO =     CelestialBody(name="Pluto",     typ="Dwarf" , mass=1/PLUTO_MASS, k_orbit=PLUTO_ORBIT)
        self.SUN =       CelestialBody(name="Sun",       typ="Star"  , mass=1/SUN_MASS, k_orbit=SUN_ORBIT)
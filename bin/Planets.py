import numpy as np
from bin.classes.AstroLi import CelestialBody, KeplerOrbit

from bin.Constants import Constants

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

class JupiterTrojan:
    """ An object containing JupiterTrojan constants

    Example:
    jt = JupiterTrojan()
    jt.HEKTOR.mass
    >> returns 624Hektors mass

    """

    def __init__(self):

        c = Constants()

        HEKTOR_ORBIT =      KeplerOrbit(5.263860462319051, 0.02262934361434583, 18.15609763764343, O=342.7878018796097, w=183.4362002247521, M=219.6678471728509)
        PATROCLUS_ORBIT =   KeplerOrbit(5.213514523425411, 0.1392756374248859 , 22.05500571757919, O=44.34827406066443, w=307.9031641166162, M=253.5524378836765)
        AGAMEMNON_ORBIT =   KeplerOrbit(5.2786680842057, 0.06677334738627223, 21.76241413631756 , O=338.005081646165, w=80.72926870121303, M=317.6174117889609)
        ACHILLES_ORBIT = KeplerOrbit(5.209257639236226, 0.1474559921319228, 10.32000402666018 , O=316.5371322563532 , w=133.3401709381824   , M=288.3217094071221   )
        MENTOR_ORBIT =   KeplerOrbit(5.149875019213133  , 0.0709249388177112, 24.63093330091959 , O=179.6064539997949, w=130.0919943058006, M=322.5148250220448)
        PARIS_ORBIT =   KeplerOrbit(5.220365495157576, 0.1278539934994891, 27.87119193219745 , O=135.9085068625309, w=150.2172062821077, M=328.0312198494113)
        DEIPHOBUS_ORBIT =   KeplerOrbit(5.129836426762841, 0.04537694142865151, 26.91437490232923, O=283.709921805651, w=359.2947329235005, M=331.9023079734244)
        ANEAS_ORBIT =   KeplerOrbit(5.221817774284521, 0.1053622430042724, 16.66200494188332, O=247.327318770212, w=50.60962526403392, M=321.4638229800369)
        DIOMEDES_ORBIT =   KeplerOrbit(5.207374114782269, 0.04414565151988255, 20.47575363906068, O=315.7779916330042, w=128.9855959977667, M=319.4954458639253)
        ODYSSEUS_ORBIT =   KeplerOrbit(5.244706312295635, 0.09099918800152354, 3.137832913162161, O=221.2416399301698, w=236.4977721243857, M=266.1390478900109 )



        HEKTOR_MASS =       7.90E18/c.M_sun #(Marchis et al, 2014)
        PATROCLUS_MASS =    1.36E18/c.M_sun #(Carry, 2012)
        AGAMEMNON_MASS=         1.36E18/c.M_sun # Just a guess??
        ACHILLES_MASS=       1.36E18/c.M_sun # Just a guess??
        MENTOR_MASS=        1.36E18/c.M_sun # Just a guess??
        PARIS_MASS=         1.36E18/c.M_sun # Just a guess??
        DEIPHOBUS_MASS= 1.36E18/c.M_sun # Just a guess??
        ANEAS_MASS=    1.36E18/c.M_sun # Just a guess??
        DIOMEDES_MASS= 1.36E18/c.M_sun # Just a guess??
        ODYSSEUS_MASS= 1.36E18/c.M_sun # Just a guess??                                 

        self.HEKTOR =       CelestialBody(name="624 Hektor",      typ="Jupiter Trojan", mass=HEKTOR_MASS,       k_orbit=HEKTOR_ORBIT)
        self.PATROCLUS =    CelestialBody(name="617 Patroclus",   typ="Jupiter Trojan", mass=PATROCLUS_MASS,    k_orbit=PATROCLUS_ORBIT)
        self.AGAMEMNON =    CelestialBody(name="911 Agamemnon",   typ="Jupiter Trojan", mass=AGAMEMNON_MASS,    k_orbit=AGAMEMNON_ORBIT)
        self.ACHILLES =     CelestialBody(name="588 Achilles",    typ="Jupiter Trojan", mass=ACHILLES_MASS,    k_orbit=ACHILLES_ORBIT)
        self.MENTOR =     CelestialBody(name="3451 Mentor",    typ="Jupiter Trojan", mass=MENTOR_MASS,    k_orbit=MENTOR_ORBIT)
        self.PARIS =     CelestialBody(name="3317 Paris",    typ="Jupiter Trojan", mass=PARIS_MASS,    k_orbit=PARIS_ORBIT)
        self.DEIPHOBUS =     CelestialBody(name="1867 Deiphobus",    typ="Jupiter Trojan", mass=DEIPHOBUS_MASS,    k_orbit=DEIPHOBUS_ORBIT)
        self.ANEAS =     CelestialBody(name="1172 Aneas",    typ="Jupiter Trojan", mass=ANEAS_MASS,    k_orbit=ANEAS_ORBIT)
        self.DIOMEDES =     CelestialBody(name="1437 Diomedes",    typ="Jupiter Trojan", mass=DIOMEDES_MASS,    k_orbit=DIOMEDES_ORBIT)
        self.ODYSSEUS =     CelestialBody(name="1143 Odysseus",    typ="Jupiter Trojan", mass=ODYSSEUS_MASS,    k_orbit=ODYSSEUS_ORBIT)

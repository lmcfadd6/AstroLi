
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from AstroLi import *

def plotOrbits(O1, O2, O3):
    """ Plots orbits from observations along with mean orbit and other planets
    Inputs:
    O1, O2, O3 [Observation Objs] - Observation objects with orbits included (run findOrbit() first!)
    """

    # If mass of the planets is 0:
    c = Constants()
    mu_sun = c.G

    # Define Keplar Orbit Objects from table
    #(a, e, i, O=None, w=None, f=None, w_tilde=None)
    Mercury = KeplerOrbit(0.38709893, 0.20563069, 7.00487, O=48.33167, w_tilde=77.45645)
    Venus = KeplerOrbit(0.72333199, 0.00677323, 3.39471, O=76.68069, w_tilde=131.53298)
    Mars = KeplerOrbit(1.52366231, 0.09341233, 1.85061, O=49.57854, w_tilde=336.04084)
    Jupiter = KeplerOrbit(5.20336301, 0.04839266, 1.30530, O=100.55615, w_tilde=14.75385)
    Saturn = KeplerOrbit(9.53707032, 0.05415060, 2.48446, O=113.71504, w_tilde=92.43194)
    Uranus = KeplerOrbit(19.19126393, 0.04716771, 0.76986, O=74.22988, w_tilde=170.96424)
    Neptune = KeplerOrbit(30.06896348, 0.00858587, 1.76917, O=131.72169, w_tilde=44.97135)
    Pluto = KeplerOrbit(39.48168677, 0.24880766, 17.14175, O=110.30347, w_tilde=224.06676)

    # https://ssd.jpl.nasa.gov/?sb_elem
    # Epoch J2000
    # Halley gives q instead of a
    # Ceres = KeplerOrbit(2.7653485, 0.07913825,  10.58682,  w=72.58981,  O=80.39320)
    # Halley = KeplerOrbit(0.58597811/(1-0.96714291), 0.96714291, 162.26269, w=111.33249, O=58.42008)


    mean = KeplerOrbit(np.mean([O1.a, O2.a, O3.a]), \
                        np.mean([O1.e, O2.e, O3.e]), \
                        np.mean([O1.i.deg, O2.i.deg, O3.i.deg]), \
                        O=np.mean([O1.O.deg, O2.O.deg, O3.O.deg]), \
                        w=np.mean([O1.w.deg, O2.w.deg, O3.w.deg]))

    # Extra information for plotting
    planets = [Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, O1, O2, O3, mean]
    colours = ["#947876", "#bf7d26", "#479ef5", "#fa0707", "#c79e0a", "#bdba04", "#02edd6", "#2200ff", "#a3986c", "#000000", "#000000", "#000000", "#fa0707"]
    names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto",  "Orbit1", "Orbit2", "Orbit3", "Mean Orbit"]
    NO_OF_PLANETS = len(planets)

    # The following animation code was taken from various stack overflow forums and from the 
    # example given on the matplotlib page
    # https://stackoverflow.com/questions/41602588/matplotlib-3d-scatter-animations

    # Helper function - Organizes all data of an orbit
    def make_planet(n, planet, less=False):
        data_x = []
        data_y = []
        data_z = []


        for f in np.linspace(0, 360, 1800):
            r, v = planet.orbit2HeliocentricState(mu_sun, f)
            data_x.append(r.x)
            data_y.append(r.y)
            data_z.append(r.z)

        data = np.array([data_x, data_y, data_z])
        return data

    # Updates a single planet
    def update(num, data, lines) :

        lines.set_data(data[0:2, num-1:num])
        lines.set_3d_properties(data[2,num-1:num])
        return lines

    # Updates all planets
    def update_all(num, data, lines):

        l = [None]*NO_OF_PLANETS

        for i in range(NO_OF_PLANETS):
            l[i] = update(num, data[i][0], lines[i][0])

        return l

    fig = plt.figure()
    ax = p3.Axes3D(fig)

    n = 100

    data = [None]*NO_OF_PLANETS
    lines = [None]*NO_OF_PLANETS

    # Generate planet lines
    for pp, p in enumerate(planets):
        data[pp] = [make_planet(n, p, less=True)]
        lines[pp] = [ax.plot(data[pp][0][0,0:1], data[pp][0][1,0:1], data[pp][0][2,0:1], \
                c=colours[pp], marker='o', label=names[pp])[0]]


    # Setthe axes properties
    ax.set_xlim3d([-5.0, 5.0])
    ax.set_xlabel('X [AU]')

    ax.set_ylim3d([-5.0, 5.0])
    ax.set_ylabel('Y [AU]')

    ax.set_zlim3d([-5.0, 5.0])
    ax.set_zlabel('Z [AU]')

    ax.set_title('Kepler Orbits')

    ax.scatter([0], [0], [0], c="y", marker='o', label="Sun")

    # Make cts line orbits
    for pp, planet in enumerate(planets):
        data_x = []
        data_y = []
        data_z = []
        for f in np.linspace(0, 2*np.pi, 1800):
            r, v = planet.orbit2HeliocentricState(mu_sun, f)
            data_x.append(r.x)
            data_y.append(r.y)
            data_z.append(r.z)

        ax.plot(data_x, data_y, data_z, c=colours[pp])

    # Creating the Animation object
    ani = animation.FuncAnimation(fig, update_all, n, fargs=(data, lines),
                                  interval=50, blit=False)
    plt.legend()
    plt.show()

def genRaDecPlot(O1, O2, O3, mag):
    """ Generates RA and Dec plots
    Inputs:
    O1, O2, O3 [Observation Objs] - Observation Objs (run findOrbit() first!)
    mag [float] - absolute magnitude of object
    """


    mu = c.G
    jd = 2452465.5
    OBL = Angle(23.439291111111, deg=True)

    ra_list = []
    dec_list = []
    jd_list = []
    mag_list = []
    sun_dist_list = []
    earth_dist_list = []
    sun_earth_list = []

    rx = []
    ry = []
    sx = []
    sy = []

    H = mag

    O = O3

    for jj in np.linspace(0, 120, 121):
        
        jy = (jd + jj)

        # get sun-earth vector at a given jd
        f = jd2f(jy, mu, Earth)
        R, _ = Earth.orbit2HeliocentricState(mu, f.rad)
        R = Vector3D(*R.xyz)

        # # get sun-obj vector at a given jd 
        f = jd2f(jy, mu, O.orbit)
        S, _ = O.orbit.orbit2HeliocentricState(mu, f.rad)
        S = Vector3D(*S.xyz)

        # earth-obj vector
        V = S - R
        V = V.rotate(-OBL, "x")


        ra, dec = cart2Radec(V)

        ra_sun, dec_sun = cart2Radec(S)

        m = H + 5*np.log10(S.mag()*V.mag()) - asteroidphasefunction(V, S)

        ra_list.append(ra)
        dec_list.append(dec)
        jd_list.append(jy)
        mag_list.append(m)

        rx.append(R.x)
        ry.append(R.y)
        sx.append(S.x)
        sy.append(S.y)
        sun_dist_list.append(S.mag())
        earth_dist_list.append(V.mag())
        sun_earth_list.append(R.mag())


    plt.subplot(2, 1, 1)
    plt.plot(jd_list, mag_list)
    plt.xlabel("Julian Day [days]")
    plt.ylabel("Apparent Magnitude")

    plt.subplot(2, 1, 2)
    plt.plot(dec_list, ra_list, )
    plt.ylabel("Right Ascension [deg]")
    plt.xlabel("Declination [deg]")

    plt.show()
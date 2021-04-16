########## Imports
from vpython import *

import numpy as np

from bin.Constants import Constants
from bin.Planets import Planets, JupiterTrojan
from bin.classes.AstroLi import CelestialBody, KeplerOrbit

###################

#GlowScript 3.0 VPython
# Stars interacting gravitationally
# Bruce Sherwood

LONG_INTEGRATION = True
OUT_FILE = "C:\\Users\\lmcfd\\Desktop\\jt_big_data.npy"
# sim = "jt"
sim = "asteroid"
    
### Browser Parameters
scene.width = 1400
scene.height = 700
# Display text below the 3D graphics:
scene.title = "Our Solar System"
scene.caption = """"""


### Import Planet Values
p = Planets()
jt = JupiterTrojan()

jd = 2451545.000000
# planets = [p.SUN, p.MERCURY, p.VENUS, p.EARTH, p.MARS, p.JUPITER, p.SATURN, p.URANUS, p.NEPTUNE, p.PLUTO]
planets = [p.SUN, p.JUPITER, p.SATURN, p.URANUS, p.NEPTUNE, p.PLUTO]
p.SUN.mass += p.MERCURY.mass + p.VENUS.mass + p.EARTH.mass + p.MARS.mass 
# Total number of spherical objects (stars, planets, etc)

if sim == "jt":
    jts = [jt.HEKTOR, jt.PATROCLUS, jt.AGAMEMNON, jt.ACHILLES, jt.MENTOR, jt.PARIS, jt.DEIPHOBUS, jt.ANEAS, jt.DIOMEDES, jt.ODYSSEUS]
elif sim == "asteroid":
    import random
    N_OBJECTS = 100
    jts = []
    for i in range(N_OBJECTS):

        a = 5.20336301
        mult = 6
        obj = KeplerOrbit(a*randrange(0.75, 1.25, 0.01), 0.04839266, 1.30530, O=100.55615, w_tilde=14.75385, L=(34.40438 + mult*i*360/N_OBJECTS)%360)

        obj_body  = CelestialBody(name="obj", typ="obj"  , mass=5.0279e-30, k_orbit=obj)
        jts.append(obj_body)

for jt in jts:
    jt.mass = 1e-16

Nstars = len(planets) + len(jts)

# Natural Units:
# Length: AU
# Time:   JY
# Mass:   M_Sun

# Constants
G = 4*np.pi**2
Msun = 1
Rsun = 0.05

# Lines length - L is a characteristic length to scale scene, axes, and simulation off of
L = 1

scene.range = 40*L
scene.forward = vec(-1,-1,-1)
###############################

####### Define 3D axes (Erased for clarity)
# xaxis = curve(color=color.gray(0.5), radius=3e8)
# xaxis.append(vec(0,0,0))
# xaxis.append(vec(L,0,0))
# yaxis = curve(color=color.gray(0.5), radius=3e8)
# yaxis.append(vec(0,0,0))
# yaxis.append(vec(0,L,0))
# zaxis = curve(color=color.gray(0.5), radius=3e8)
# zaxis.append(vec(0,0,0))
# zaxis.append(vec(0,0,L))
data = []
Stars = []
star_colors = [color.red, color.green, color.blue,
              color.yellow, color.cyan, color.magenta]


### Input planet parameters
mas_list = []
pos_list = []
vel_list = []
# I had to make the outer planets larger so they were visible, here's a scaling factor
if sim == "jt":
    sizes = [6]*len(jts)
elif sim == "asteroid":
    sizes = [1]*len(jts)

rad_list = [20, 12, 12, 12, 12, 12] + sizes

for planet in planets:
    m, r, v = planet.status(jd)
    mas_list.append(m)
    pos_list.append(vec(r.x, r.y, r.z))
    vel_list.append(vec(v.x, v.y, v.z))

m, r_this, v_this = p.JUPITER.status(jd)

for planet in jts:
    m, r, v = planet.status(jd)
    mas_list.append(m)
    pos_list.append(vec(r.x, r.y, r.z))
    vel_list.append(vec(v.x, v.y, v.z))
########################################

### This loop initializes the stars, this is where the initial parameters from the actual solar system go
psum = vec(0,0,0)
for i in range(Nstars):
    # Make sphere with initial position
    # Retain - how long to keep tail
    # trail_radius, size of the trail
    if i < len(planets):
        # if sim == "jt":
        #     star = sphere(pos=pos_list[i], make_trail=True, retain=1500, trail_radius=0.001*rad_list[i])
        # elif sim == "asteroid":
        star = sphere(pos=pos_list[i], make_trail=False)

    else:
        star = sphere(pos=pos_list[i], make_trail=False)
    R = Rsun
    star.radius = R*rad_list[i]

    star.mass = mas_list[i]
    star.momentum = mas_list[i]*vel_list[i]

    # Assign colors from list (loop over)
    star.color = star.trail_color = star_colors[i % 6]
    Stars.append( star )
    # Total momentum of the system
    psum = psum + star.momentum
    if i == 0: continue

    # data.append([i, star.pos, star.momentum/star.mass, G*(1 + star.mass)])

#make total initial momentum equal zero

# Ignore this, my momenta are fine
# for i in range(Nstars):
#     Stars[i].momentum = Stars[i].momentum - psum/Nstars

dt = (11)*0.01 # <- Thanks Cole! So that 1% of Jupiter orbit is timestep
hitlist = []

def forceOnStar(si, N, i, r):
    
    F = vec(0,0,0)
    pos1 = r
    m1 = si.mass
    radius = si.radius

    M_tot = Msun
    F_indirect = vec(0,0,0)
    for j in range(N):
        if j == 0: continue
        sj = Stars[j]
        dis = sj.pos

        F_indirect = F_indirect + (G*sj.mass/(mag2(dis)**1.5))*(dis)
        M_tot += sj.mass

    # this is the sum of all forces with i, j with the case i != j
    for j in range(N):
        if i == j: continue
        sj = Stars[j]

        r = sj.pos - pos1
        rmag2 = mag2(r)
        if rmag2 <= (radius+sj.radius)**2: hitlist.append([i,j])
        # Calculate force (rmag2 is r^2, rmag2^1.5 is r^3, so we get r/r^3)
        # Standard Newton force
        try:
            F = F + (G*m1*sj.mass/(rmag2**1.5))*r
        except ZeroDivisionError:
            # cheat and use machine eps if it would be 0, so it at least doesn't break
            F = F + (G*m1*sj.mass/(7/3-4/3-1))*r

    F = F - F_indirect/M_tot*m1

    return F


def computeForces():
    global hitlist, Stars
    hitlist = []
    N = len(Stars)
    # Calculate the force on every star (find new momentum so that we change the position next time step)
    for i in range(N):

        si = Stars[i]
        # don't integrate the sun!
        if i == 0:
            si.pos = vec(0, 0, 0)
            si.momentum = vec(0, 0, 0)
            continue

        

        # # Change momentum vector - Euler's integration here too
        if LONG_INTEGRATION:

            p_1 = si.momentum
            r = si.pos

            F = forceOnStar(si, N, i, r) # F(t_0, p_0)

            k_1 = F

            p_2 = p_1 + F*dt/2 
            F = forceOnStar(si, N, i, r)

            k_2 = F

            r = r + p_2*(dt/2/si.mass)
            F = forceOnStar(si, N, i, r)

            k_3 = F

            p_3 = p_2 + F*dt/2
            r = r + p_3*(dt/2/si.mass)

            F = forceOnStar(si, N, i, r)

            k_4 = F

            si.momentum = si.momentum + dt/6*(k_1 + 2*k_2 + 2*k_3 + k_4)
            si.pos = si.pos + dt/6*(p_1 + 2*p_2 + 2*p_3 + si.momentum)/si.mass

        else:

            F = forceOnStar(si, N, i, si.pos)

            si.momentum = si.momentum + F*dt

            si.pos = si.pos + si.momentum*(dt/si.mass)


        # planet number, pos, vel, mu
        # data.append([i, r_this, v_this, G*(Stars[0].mass + si.mass)])
        v = si.momentum/si.mass

        # Did not work nicely with vectors because of some pickling thing - so instead we have this mess!
        data.append([i, si.pos.x, si.pos.y, si.pos.z, v.x, v.y, v.z, G*(1 + si.mass)])

counter = 0
# Main animation loop
while True:

    counter += 1
    # Rate of viewing (not timestep)
    rate(1000)
    
    # Compute all forces on all stars
    computeForces()

    if counter%10000 == 0:
        np.save(OUT_FILE, np.array(data))
    # # Having updated all momenta, now update all positions
    # for star in Stars:
    #     # don't calculate collided star twice
    #     if star is None: continue
        
    #     # Euler method of integration
        

    # # If any collisions took place, merge those stars
    # hit = len(hitlist)-1

    ### My code would get stuck in this loop for some reason, since we don't need collision, I just omitted this section
    # # Collisions take place here
    # while hit > 0:
    #     s1 = Stars[hitlist[hit][0]]
    #     s2 = Stars[hitlist[hit][1]]
    #     if (s1 is None) or (s2 is None): continue
    
    #     mass = s1.mass + s2.mass
    #     momentum = s1.momentum + s2.momentum
    #     pos = (s1.mass*s1.pos + s2.mass*s2.pos) / mass

    #     # Cool way of merging colours based off of weighted colour
    #     s1.color = s1.trail_color = (s1.mass*s1.color + s2.mass*s2.color) / mass
    #     R = Rsun
        
    #     # Remove collided star
    #     s1.clear_trail()
    #     s2.clear_trail()
    #     s2.visible = False
        
    #     s1.mass = mass
    #     s1.momentum = momentum
    #     s1.pos = pos
    #     s1.radius = R
    #     Stars[hitlist[hit][1]] = None
    #     hit -= 1

from vpython import *

import numpy as np

from bin.Constants import Constants
from bin.Planets import Planets

#GlowScript 3.0 VPython

# Stars interacting gravitationally

# Bruce Sherwood

## Screen size in browser
scene.width = scene.height = 600

# Display text below the 3D graphics:
scene.title = "Stars interacting gravitationally"

scene.caption = """To rotate "camera", drag with right button or Ctrl-drag.
To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
  On a two-button mouse, middle is left + right.
To pan left/right and up/down, Shift-drag.
Touch screen: pinch/extend to zoom, swipe or two-finger rotate.

Replot the web page to restart with different (random) initial conditions."""

p = Planets()

jd = 2451545.000000
planets = [p.SUN, p.MERCURY, p.VENUS, p.EARTH, p.MARS, p.JUPITER, p.SATURN, p.URANUS, p.NEPTUNE, p.PLUTO]
# Total number of spherical objects (stars, planets, etc)
Nstars = len(planets)

# Natural Units:
# Length: AU
# Time:   JD
# Mass:   M_Sun

# Constants
G = 4*np.pi**2
Msun = 1
Rsun = 0.05

# Lines length
L = 1

scene.range = 2*L
scene.forward = vec(-1,-1,-1)

# # Define 3D axes
# xaxis = curve(color=color.gray(0.5), radius=3e8)
# xaxis.append(vec(0,0,0))
# xaxis.append(vec(L,0,0))
# yaxis = curve(color=color.gray(0.5), radius=3e8)
# yaxis.append(vec(0,0,0))
# yaxis.append(vec(0,L,0))
# zaxis = curve(color=color.gray(0.5), radius=3e8)
# zaxis.append(vec(0,0,0))
# zaxis.append(vec(0,0,L))

Stars = []
star_colors = [color.red, color.green, color.blue,
              color.yellow, color.cyan, color.magenta]



mas_list = []
pos_list = []
vel_list = []
rad_list = [1, 1, 1, 1, 1, 2, 3, 4, 5, 6]

for planet in planets:
    m, r, v = planet.status(jd)
    mas_list.append(m)
    pos_list.append(vec(r.x, r.y, r.z))
    vel_list.append(vec(v.x, v.y, v.z))


### This loop initializes the stars, this is where the parameters from the actual solar system go
psum = vec(0,0,0)
for i in range(Nstars):
    # Make sphere with initial position
    star = sphere(pos=pos_list[i], make_trail=True, retain=1500, trail_radius=0.01*rad_list[i])
    R = Rsun
    star.radius = R*rad_list[i]
    # Program guesses mass from size (this will be changed with actual mass)
    star.mass = mas_list[i]

    # Random momentum assinged (this will also be changed)
    star.momentum = mas_list[i]*vel_list[i]
    star.color = star.trail_color = star_colors[i % 6]
    Stars.append( star )
    # Total momentum of the system
    psum = psum + star.momentum

#make total initial momentum equal zero

# Ignore this, my momenta are fine
# for i in range(Nstars):
#     Stars[i].momentum = Stars[i].momentum - psum/Nstars

dt = (88/365)*0.01*0.01 # <- Thanks Cole!
hitlist = []

def computeForces():
    global hitlist, Stars
    hitlist = []
    N = len(Stars)
    # Calculate the force on every star (find new momentum so that we change the position next time step)
    for i in range(N):
        si = Stars[i] 
        if si is None: continue
        F = vec(0,0,0)
        pos1 = si.pos
        m1 = si.mass
        radius = si.radius

        # this is the sum of all forces with i, j with the case i != j
        for j in range(N):
            if i == j: continue
            sj = Stars[j]
            if sj is None: continue
            r = sj.pos - pos1
            rmag2 = mag2(r)
            if rmag2 <= (radius+sj.radius)**2: hitlist.append([i,j])
            # Calculate force (rmag2 is r^2, rmag2^1.5 is r^3, so we get r/r^3)
            F = F + (G*m1*sj.mass/(rmag2**1.5))*r

        # Change momentum vector
        si.momentum = si.momentum + F*dt


# Main animation loop
while True:

    # Rate of viewing (not timestep)
    rate(100000)
    
    # Compute all forces on all stars
    computeForces()

    # Having updated all momenta, now update all positions
    for star in Stars:
        # don't calculate collided star twice
        if star is None: continue
        
        # Euler method of integration
        star.pos = star.pos + star.momentum*(dt/star.mass)

    # If any collisions took place, merge those stars
    hit = len(hitlist)-1

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

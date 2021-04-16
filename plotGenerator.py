import numpy as np
import matplotlib.pyplot as plt

from bin.Helper import orbitalElements
from bin.classes.AstroLi import *

IN_FILE = "C:\\Users\\lmcfd\\Desktop\\jt_big_data.npy"

data = np.load(IN_FILE, allow_pickle=True)
SKIP = 10

dt = (11)*0.01

jupiter_data = []
saturn_data = []
uranus_data = []
neptune_data = []
pluto_data = []
jt_data = []
r_jupiter = []
r_hektor = []


for line in data:

	obj = line[0]


	# Ignore Sun
	if obj == 0: continue


	r = Vector3D(line[1], line[2], line[3])
	v = Vector3D(line[4], line[5], line[6])
	mu = line[7]

	a, e, i, o, f, w = orbitalElements(mu, r, v, write=False)

	try:
		o = o%360
	except TypeError:
		o = None
	try:
		w = w%360
	except TypeError:
		w = None

	if obj == 1:
		jupiter_data.append([a, e, i, o, w, f])
		r_jupiter.append(r)
	elif obj == 2:
		saturn_data.append([a, e, i, o, w, f])
	elif obj == 3:
		uranus_data.append([a, e, i, o, w, f])
	elif obj == 4:
		neptune_data.append([a, e, i, o, w, f])
	elif obj == 5:
		pluto_data.append([a, e, i, o, w, f])
	else:
		jt_data.append([a, e, i, o, w, f])
		r_hektor.append(r)

dt = np.arange(0, len(jupiter_data), 1)*dt


jupiter_data = np.array(jupiter_data)
saturn_data = np.array(saturn_data)
neptune_data = np.array(neptune_data)
uranus_data = np.array(uranus_data)
pluto_data = np.array(pluto_data)
jt_data = np.array(jt_data)

labels = ["Semimajor Axis, a [AU]", "Eccentricity, e", "Inclination, i, [deg]", "Longitude of Acending Node, Omega, [deg]", \
			"Argument of Periapsis, omega, [deg]"]

# for i in range(5):
# 	plt.subplot(3, 2, i+1)

# 	plt.plot(t, jupiter_data[:, i], label="Jupiter")
# 	plt.plot(t, saturn_data[:, i], label="Saturn")
# 	plt.plot(t, uranus_data[:, i], label="Uranus")
# 	plt.plot(t, neptune_data[:, i], label="Neptune")
# 	# plt.plot(t[::SKIP], pluto_data[::SKIP, i], label="Pluto")
# 	plt.xlabel("Time [yrs]")
# 	plt.ylabel(labels[i])
# 	plt.legend()

# plt.show()

t = np.repeat(np.arange(0, len(jupiter_data), 1)*dt, obj-5)

# for i in range(5):
# 	plt.subplot(3, 2, i+1)
# 	plt.scatter(t, jt_data[:, i], label="Trojans")
# 	# plt.plot(dt, jupiter_data[:, i], label="Jupiter")
# 	plt.xlabel("Time [yrs]")
# 	plt.ylabel(labels[i])
# 	plt.legend()

# plt.show()

# plt.plot(t, r_jupiter, label="Jupiter")
# plt.plot(t, r_hektor, label="Hektor")
# plt.xlabel("Time [yrs]")
# plt.ylabel("Distance from the Sun [AU]")
# plt.legend()
# plt.show()
dif = []
# for rr in range(len(r_jupiter)):
az_j = np.arctan2(r_jupiter[-1].x, r_jupiter[-1].y)
for i in range(100):
	az_h = np.arctan2(r_hektor[-1 - i].x, r_hektor[-1 - i].y)

	dif.append(np.degrees(az_j - az_h)%360)

plt.hist(dif, bins=18)
plt.xlabel("Angular Separation [deg]")
plt.ylabel("Number of Asteroids/100")
plt.show()
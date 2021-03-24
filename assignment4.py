from bin.Constants import Constants
from bin.Planets import Planets

p = Planets()

if __name__ == '__main__':

	jd = 2451545.000000

	p.MERCURY.status(jd, write=True)
	p.VENUS.status(jd, write=True)
	p.EARTH.status(jd, write=True)
	p.MARS.status(jd, write=True)
	p.JUPITER.status(jd, write=True)
	p.SATURN.status(jd, write=True)
	p.URANUS.status(jd, write=True)
	p.NEPTUNE.status(jd, write=True)
	p.PLUTO.status(jd, write=True)
	
	jd = 2452545.000000

	p.MERCURY.status(jd, write=True)
	p.VENUS.status(jd, write=True)
	p.EARTH.status(jd, write=True)
	p.MARS.status(jd, write=True)
	p.JUPITER.status(jd, write=True)
	p.SATURN.status(jd, write=True)
	p.URANUS.status(jd, write=True)
	p.NEPTUNE.status(jd, write=True)
	p.PLUTO.status(jd, write=True)
	
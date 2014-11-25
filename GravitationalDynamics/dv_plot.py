import numpy
from matplotlib import pyplot

from amuse.units import units, nbody_system
from amuse.plot import scatter, xlabel, ylabel

converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.parsec)
G_si = converter.to_si(nbody_system.G)
milkyway_mass = 1.26e12 | units.MSun
M = milkyway_mass

def escape_velocity(impact_parameter):
    return ((2 * G_si * M)/impact_parameter).sqrt()

def v_infinity(escape_v, v):
    return (v**2 - escape_v**2).sqrt()



pyplot.title("d,v-plane")
for i in numpy.arange(0.01,1,0.01):
    d = 100000|units.parsec * i
    esc_velocity = escape_velocity(d)
    for j in numpy.arange(0.01,1,0.01):
        v = 10000*j|units.parsec/units.Myr
        if v < esc_velocity:
            scatter(d,v, color='blue', s=2)
        else:
            scatter(d,v, color='red', s=2)

xlabel("d")
ylabel("v")
pyplot.xlim(0,10000000)
pyplot.ylim(0,10000)
pyplot.savefig("plots/dv_plot")            



import numpy
from matplotlib import pyplot
from math import pi
import pickle 

from amuse.community.bhtree.interface import BHTree
from amuse.lab import new_plummer_model
from amuse.units import units, nbody_system
from amuse.datamodel import Particle
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import zero

from find_bound_mass import find_bound_bodies, find_bound_bodies_using_energies

R = 6000 | units.parsec

def plot_cluster(stars, filename, t, rcl):
    filename += ".png"
    #lim = 10*stars.center_of_mass().length().value_in(stars.x.unit)
    lim = (R * 2).value_in(units.m)
    m = 1 + 3.0*stars.mass/min(stars.mass)

    pyplot.title("Cluster at distance "+str(R)+" with cluster radius "+
            str(rcl)+"\nat t="+str(t))

    scatter(stars.x, stars.y)  # s size (in point^2)
    xlabel("X")
    ylabel("Y")
    pyplot.xlim(-lim, lim)
    pyplot.ylim(-lim, lim)
    pyplot.savefig(filename)
    pyplot.clf()


def nbody_integrator(Ncl, mcl, rcl, t_end, n_steps, algorithm=BHTree, 
        timestep=None, blackhole=False):
    converter = nbody_system.nbody_to_si(mcl, rcl)
    bodies = new_plummer_model(Ncl, convert_nbody=converter)

    #estimate of milky way mass by "Mass models of the Milky Way", McMillan
    blackhole_mass = 1.26e12 | units.MSun
    blackhole = Particle(mass=blackhole_mass)
    blackhole.position = [0,0,0] | units.m

    cluster_velocity = [0,0,0] | units.m / units.s
    cluster_position = [0,0,0] | units.parsec
    cluster_position[0] = R
    G_si = converter.to_si(nbody_system.G)
    V = (G_si * blackhole_mass/R).sqrt().as_quantity_in(units.m/units.s)
    cluster_velocity[1] = V
    bodies.move_to_center()
    bodies.velocity += cluster_velocity
    bodies.position += cluster_position
    bodies.add_particle(blackhole)
        

    gravity = BHTree(converter)
    gravity.particles.add_particles(bodies)
    channel_from_gravity_to_framework = gravity.particles.\
        new_channel_to(bodies)

    time = zero
    dt = t_end / float(n_steps)
    x = 0
    while time < t_end:
        plot_cluster(bodies, "blackhole_plots/"+
                str(radius.value_in(units.parsec))+"_"+str(x),time, rcl)
        time += dt
        gravity.evolve_model(time)
        channel_from_gravity_to_framework.copy()
        x+=1
    plot_cluster(bodies, "blackhole_plots/"+str(x),time, rcl)
    gravity.stop()
    return bodies


def find_cluster_bound_mass_fraction(bodies, energy=False):
    cluster_mass = 10**7 | units.MSun
    cluster = bodies[:-1]
    if energy:
        bound_bodies = find_bound_bodies_using_energies(cluster)
    else:
        bound_bodies = find_bound_bodies(cluster)
    bound_mass = bound_bodies.total_mass()
    return bound_mass / cluster.total_mass()





if __name__ in '__main__':
    options = {}
    options['mcl'] = 10**7 | units.MSun
    options['n_steps'] = 200
    options['algorithm'] = BHTree
    N = 1024
    options['Ncl'] = N
    

    converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.parsec)
    G_si = converter.to_si(nbody_system.G)
    blackhole_mass = 1.26e12 | units.MSun
    t_end = 2* pi * (R**3/(G_si*blackhole_mass)).sqrt()
    options['t_end'] = t_end *2

    datapoints = []

    for i in xrange(10, 101,10):
        radius = i | units.parsec
        options['rcl'] = radius
        bodies = nbody_integrator(**options)
        hop_mass_fraction = find_cluster_bound_mass_fraction(bodies)
        energy_mass_fraction = find_cluster_bound_mass_fraction(bodies, energy=True)

        data={}
        data["distance"] = R
        data["radius"] = radius
        data["hop_mass_fraction"] = hop_mass_fraction
        data["energy_mass_fraction"] = energy_mass_fraction
        print radius, hop_mass_fraction, energy_mass_fraction
        datapoints.append(data)
    pickle.dump(datapoints, open("bound_mass.dat", "wb"))




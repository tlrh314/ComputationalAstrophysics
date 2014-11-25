import numpy
from matplotlib import pyplot
from math import pi
import pickle 

from amuse.community.bhtree.interface import BHTree
from amuse.lab import new_plummer_model
from amuse.units import units, nbody_system
from amuse.datamodel import Particle
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import zero, Kepler

from find_bound_mass import find_bound_bodies, find_bound_bodies_using_energies


def plot_cluster(stars, filename, t, rcl, V):
    filename += ".png"
    #lim = 10*stars.center_of_mass().length().value_in(stars.x.unit)
    lim = (R * 2).value_in(units.m)
    m = 1 + 3.0*stars.mass/min(stars.mass)

    pyplot.title("Cluster at distance "+str(R)+" with cluster radius "+
            str(rcl)+"\nvelocity "+str(V)+"\nat t="+str(t))

    scatter(stars.x, stars.y)  # s size (in point^2)
    xlabel("X")
    ylabel("Y")
    pyplot.xlim(-lim, lim)
    pyplot.ylim(-lim, lim)
    pyplot.savefig(filename)
    pyplot.clf()


def nbody_integrator(Ncl, mcl, rcl, t_end, n_steps, escape_velocity_fraction, R):
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
    escape_v = (2*G_si*blackhole_mass/R).sqrt().as_quantity_in(units.m/units.s)
    V = escape_v * escape_velocity_fraction
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
    base_path = "encounter_plots/"+str(R.value_in(units.parsec))+"_"+\
        str(escape_velocity_fraction)+"_"
    while time < t_end:
        plot_cluster(bodies, base_path+str(x),time, rcl, V)
        time += dt
        gravity.evolve_model(time)
        channel_from_gravity_to_framework.copy()
        x+=1
    plot_cluster(bodies, base_path+str(x),time, rcl, V)

    gravity.stop()
    return V, bodies


def find_cluster_bound_mass_fraction(bodies, energy=False):
    cluster_mass = 10**7 | units.MSun
    cluster = bodies[:-1]
    if energy:
        bound_bodies = find_bound_bodies_using_energies(cluster)
    else:
        bound_bodies = find_bound_bodies(cluster)
    bound_mass = bound_bodies.total_mass()
    return bound_mass / cluster.total_mass()





def find_kepler_parameters(bodies, converter):
    kep = Kepler(converter)
    kep.initialize_code()
    bh = bodies[-1]
    cluster = bodies[:-1]
    pos = bh.position - cluster.center_of_mass()
    vel = bh.velocity - cluster.center_of_mass_velocity()
    kep.initialize_from_dyn(bh.mass + cluster.total_mass(), pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    kep.stop()
    return a,e


if __name__ in '__main__':
    options = {}
    options['mcl'] = 10**7 | units.MSun
    options['n_steps'] = 200
    N = 1024
    options['Ncl'] = N
    

    converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.parsec)
    G_si = converter.to_si(nbody_system.G)
    blackhole_mass = 1.26e12 | units.MSun

    datapoints = []

    for i in [3000, 6000, 12000, 18000]:
        R = i | units.parsec
        t_end = 2* pi * (R**3/(G_si*blackhole_mass)).sqrt()
        options['t_end'] = t_end *2
        options['R'] = R
        radius = 0.005 * R
        options['rcl'] = radius
        for j in [0.8, 1.0, 1.2]:
            options["escape_velocity_fraction"] = j
            V, bodies = nbody_integrator(**options)
            hop_mass_fraction = find_cluster_bound_mass_fraction(bodies)
            a,e = find_kepler_parameters(bodies, converter)

            data={}
            data["distance"] = R
            data["radius"] = radius
            data["hop_mass_fraction"] = hop_mass_fraction
            data["escape_velocity_fraction"] = j
            data["velocity"] = V
            data["kepler_params"] = (a,e)
            print radius, hop_mass_fraction
            print a
            print e
            datapoints.append(data)
    pickle.dump(datapoints, open("bound_mass_encounter.dat", "wb"))




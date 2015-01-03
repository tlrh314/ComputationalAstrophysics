import numpy

from amuse.units import units, constants, nbody_system
from amuse.units.optparse import OptionParser

from amuse.support.console import set_printing_strategy
from amuse.io import write_set_to_file
from amuse.datamodel.particles import Channels

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ic.plummer import new_plummer_model

from amuse.couple.bridge import Bridge
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree

def create_cluster(Ncl, Mcl, Rcl):
    converter = nbody_system.nbody_to_si(Mcl, Rcl)
    particles = new_plummer_model(Ncl, convert_nbody=converter)
    particles.move_to_center()
    return particles, converter

def create_giant_molecular_clouds(Ngas, Mgas, Rgas, Nclouds):
    converter = nbody_system.nbody_to_si(Mgas, Rgas)
    clouds = []
    for i in xrange(Nclouds):
        particles = molecular_cloud(targetN=Ngas, convert_nbody=converter).result
        particles.u = 1. | units.ms**2
        particles.velocity *= 0.6
        particles.move_to_center()
        clouds.append(particles)
    return clouds, converter

def put_in_orbit(cluster, clouds, d, v_inf):
    gas = clouds[0]
    r = 2 * (abs(cluster.position).max() + abs(gas.position).max())
    r = max(r, 2*d)

    mu = constants.G * (cluster.mass.sum() + gas.mass.sum())
    a = -mu/v_inf**2
    e = -d/a + 1.
    v = (mu * (2./r - 1./a)).sqrt()

    true_anomaly = numpy.arccos(a*(1.-e**2)/(e*r) - 1./e)
    flight_path_angle = numpy.arctan(e*numpy.sin(true_anomaly)/
            (1 + e*numpy.cos(true_anomaly)))
    approach_angle = numpy.pi/2. - flight_path_angle

    cluster.x += r
    cluster.vx -= v * numpy.cos(approach_angle)
    cluster.vy += v * numpy.sin(approach_angle)

    asymptote_intersection_distance = r * numpy.tan(approach_angle)

    # angle between approach and departure (along asymptotes) is 2*theta
    theta = numpy.arccos(1./e)
    outgoing_angle = 2 * theta - numpy.pi/2. - approach_angle


    # place additional clouds on a direct path of the outgoing cluster 
    for i, cloud in enumerate(clouds[1:]):
        cloud.x +=  numpy.sin(outgoing_angle) * 5 * r * (i+1)
        cloud.y +=  asymptote_intersection_distance + \
                numpy.cos(outgoing_angle) * 5 * r * (i+1)
        gas.add_particles(cloud)


    t_end = (2 * r / v) * len(clouds) + (len(clouds)-1) * (5*r)/v
    return cluster, gas, t_end

def setup_codes(cluster, gas, nbody_converter, gas_converter):
    gravity = BHTree(nbody_converter)
    gravity.particles.add_particles(cluster)

    hydro = Fi(gas_converter)
    hydro.parameters.eps_is_h_flag = True
    hydro.parameters.timestep = 0.005 | units.Myr
    hydro.particles.add_particles(gas)

    bridge = Bridge(timestep=0.01|units.Myr, use_threading=False)
    bridge.add_system(gravity, (hydro,))
    bridge.add_system(hydro, (gravity,))
    
    channels = Channels()
    channels.add_channel(gravity.particles.new_channel_to(cluster))
    channels.add_channel(hydro.particles.new_channel_to(gas))

    return gravity, hydro, bridge, channels

def evolve_cluster_and_cloud(Ncl, Mcl, Rcl, Ngas, Mgas, Rgas,
        d, v_inf, dt, filename, Nclouds):
    cluster, nbody_converter = create_cluster(Ncl, Mcl, Rcl)
    clouds, gas_converter = create_giant_molecular_clouds(Ngas, Mgas, 
            Rgas, Nclouds)
    cluster, gas, t_end = put_in_orbit(cluster, clouds, d, v_inf)
    print "t end", t_end

    gravity, hydro, bridge, channels = setup_codes(cluster, gas,
            nbody_converter, gas_converter)

    time = 0. | units.Myr
    write_set_to_file((cluster.savepoint(time), gas.savepoint(time)), filename,
            'amuse', names=("cluster", "gas"), append_to_file=False)

    while time < (t_end - dt/2.):
        time += dt
        print "evolving to", time
        bridge.evolve_model(time)
        channels.copy()
        write_set_to_file((cluster.savepoint(time), gas.savepoint(time)),
                filename, 'amuse', names=("cluster", "gas"))

if __name__ == "__main__":
    set_printing_strategy("custom",
            preferred_units = [units.MSun, units.parsec, units.Myr, units.kms])

    cluster_stars = 128
    cluster_mass = 10**6|units.MSun 
    cluster_radius = 10|units.parsec
    gas_particles = 1024
    gas_mass = 10**7|units.MSun
    gas_radius = 100|units.parsec
    impact_parameter = 150|units.parsec
    v_infinity = 1|units.kms
    timestep = 0.2|units.Myr
    number_of_clouds = 2
    filename = "test.amuse"
    evolve_cluster_and_cloud(cluster_stars, cluster_mass, cluster_radius, 
            gas_particles, gas_mass, gas_radius, impact_parameter, 
            v_infinity, timestep, filename, number_of_clouds)





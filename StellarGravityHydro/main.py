from amuse.lab import *
from amuse.units.optparse import OptionParser
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.couple.bridge import Bridge
from amuse.datamodel.particles import Channels


import numpy

def create_cluster(Ncl, Mcl, Rcl):
    converter = nbody_system.nbody_to_si(Mcl, Rcl)
    particles = new_plummer_model(Ncl, convert_nbody=converter)
    particles.move_to_center()
    return particles, converter

def create_giant_molecular_cloud(Ngas, Mgas, Rgas):
    converter = nbody_system.nbody_to_si(Mgas, Rgas)
    particles = molecular_cloud(targetN=Ngas, convert_nbody=converter).result
    particles.u = 1. | units.ms**2
    particles.velocity *= 0.6
    particles.move_to_center()
    return particles, converter

def put_in_orbit(cluster, gas, d, v_inf):
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

    t_end = 2 * r / v
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

    #####################
    stellar = SSE()
    stars = stellar.particles.add_particles(cluster)

    from_stellar_evolution_to_model \
        = stellar.particles.new_channel_to(cluster)
    from_stellar_evolution_to_model.copy_attributes(["mass"])
    #####################

    return gravity, hydro, bridge, channels, stellar, from_stellar_evolution_to_model


def evolve_cluster_and_cloud(Ncl, Mcl, Rcl, Ngas, Mgas, Rgas,
        d, v_inf, dt, filename):
    cluster, nbody_converter = create_cluster(Ncl, Mcl, Rcl)
    gas, gas_converter = create_giant_molecular_cloud(Ngas, Mgas, Rgas)
    cluster, gas, t_end = put_in_orbit(cluster, gas, d, v_inf)
    print "t end", t_end

    gravity, hydro, bridge, channels, stellar, stellar_to_cluster = setup_codes(cluster, gas,
            nbody_converter, gas_converter)

    time = 0. | units.Myr
    write_set_to_file((cluster.savepoint(time), gas.savepoint(time)), filename,
            'amuse', names=("cluster", "gas"), append_to_file=False)

    while time < (t_end - dt/2.):
        time += dt
        print "evolving to", time
        stellar.evolve_model()
        stellar_to_cluster.copy_attributes(["mass", "radius"])
        bridge.evolve_model(time)
        channels.copy()
        write_set_to_file((cluster.savepoint(time), gas.savepoint(time)),
                filename, 'amuse', names=("cluster", "gas"))

def parse_options():
    parser = OptionParser()
    parser.add_option("-N", dest="Ncl", type="int",
            default = 128, help="number of stars [%default]")
    parser.add_option("-M", dest="Mcl", unit=units.MSun,
            type="float", default = 10**6|units.MSun,
            help="clustermass [%default]")
    parser.add_option("-R", dest="Rcl", unit=units.parsec,
            type="float", default = 10|units.parsec,
            help="cluster half-mass radius [%default]")
    parser.add_option("-n", dest="Ngas", type="int",
            default = 1024, help="number of gas particles [%default]")
    parser.add_option("-m", dest="Mgas", unit=units.MSun,
            type="float", default = 10**7|units.MSun,
            help="gas mass [%default]")
    parser.add_option("-r", dest="Rgas", unit=units.parsec,
            type="float", default = 100|units.parsec,
            help="gas radius [%default]")
    parser.add_option("-d", dest="d", unit=units.parsec,
            type="float", default = 150|units.parsec,
            help="Impact parameter of the encounter [%default]")
    parser.add_option("-v", dest="v_inf", unit=units.kms,
            type="float", default = 1|units.kms,
            help="Cluster velocity at infinity [%default]")
    parser.add_option("-T", dest="dt", unit=units.Myr,
            type="float", default = 0.2|units.Myr,
            help="Snapshot output timestep [%default]")
    parser.add_option("-f", dest="filename",
            type="string", default = "test.amuse",
            help="The file in which to save the snapshots [%default]")

    options, arguments = parser.parse_args()
    return options.__dict__

if __name__ == "__main__":
    set_printing_strategy("custom",
            preferred_units = [units.MSun, units.parsec, units.Myr, units.kms])

    options = parse_options()
    evolve_cluster_and_cloud(**options)

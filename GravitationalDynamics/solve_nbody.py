#!/usr/bin/env amuse
"""
solve_nbody.py

Computational Astrophysics / Gravitational Dynamics
Simulating Ultra Compact Dwarf Galaxies with AMUSE

Timo Halbesma, s1603221.
Shebaz Sultan, s1617451.

November 25, 2014. Version 1.0

"""

import traceback

from numpy import isnan

# from amuse.lab import *
from amuse.community.bhtree.interface import BHTree
from amuse.lab import new_plummer_model, write_set_to_file
from amuse.lab import set_printing_strategy, zero
from amuse.units import units, nbody_system


def nbody_integrator(Ncl, mcl, rcl, t_end, n_steps, algorithm=BHTree):
    converter = nbody_system.nbody_to_si(mcl, rcl)
    bodies = new_plummer_model(Ncl, convert_nbody=converter)

    verbose = True # For AMUSE functions this is called redirection.
    algorithm_name = str(algorithm).split('.')[-1][:-2]
    runtime_of_N_and_t = numpy.zeros((len(stars), len(t_end)),
                                     dtype=numpy.float64)
    dE_of_N_and_t = numpy.zeros((len(stars), len(t_end)),
                                dtype=numpy.float64)

    # Allow selecting Stellar Dynamics code trough a function argument.
    try:
        if algorithm_name == 'Hermite':
            gravity = algorithm(converter, number_of_workers=8)
        else:
            gravity = algorithm(converter)
    except Exception, e:  # Too generic, but I don't know what Error to expect.
        # raise Exception
        print "Failure for {0}".format(algorithm_name)
        print Exception, e
        print traceback.format_exc()
        return None
    else:
        print "\nUsing algorithm '{0}', Ncl={1}, t_end={2}"\
            .format(algorithm_name, Ncl, t_end)


    # Assignment 1B
    # gravity.parameters.timestep = 0.03125

    gravity.particles.add_particles(bodies)
    channel_from_gravity_to_framework = gravity.particles.\
        new_channel_to(bodies)

    write_set_to_file(bodies.savepoint(0.0 | t_end.unit),
                      "data/nbody_{0}_{1}_{2}.hdf5"
                      .format(algorithm_name, Ncl, t_end),
                      "hdf5", append_to_file=False)

    Etot_init = gravity.kinetic_energy + gravity.\
        potential_energy

    time = zero
    dt = t_end / float(n_steps)
    while time < t_end:
        Etot_prev = Etot_init
        time += dt

        gravity.evolve_model(time)
        channel_from_gravity_to_framework.copy()
        write_set_to_file(bodies.savepoint(time),
                          "data/nbody_{0}_{1}_{2}.hdf5"
                          .format(algorithm_name, Ncl, t_end),
                          "hdf5")

        Ekin = gravity.kinetic_energy
        Epot = gravity.potential_energy
        Etot = Ekin + Epot

        # I had some concerns writing to stdout could lower the runtime
        # (wall-clock) because it tends to be slow.
        if verbose:
            print "T =", time, "M =", bodies.mass.sum(), "E =", Etot, "Q =",\
                Ekin / Epot,
            print "dE =", (Etot_init - Etot) / Etot, "ddE =",\
                (Etot_prev - Etot) / Etot

    gravity.stop()

    dE = (Etot_init - Etot) / Etot
    if isnan(dE):  # numpy function
        return 0
    else:
        return dE


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="Ncl", type="int", default=100,
                      help="number of stars [%default]")
    result.add_option("-t", unit=units.Myr, dest="t_end", type="float",
                      default=1 | units.Myr,
                      help="end time of the simulation [%default]")
    result.add_option("-n", dest="n_steps", type="float", default=100,
                      help="number of output steps [%default]")
    result.add_option("-m", unit=units.parsec, dest="mcl", type="float",
                      default=10**7 | units.MSun,
                      help="cluster mass [%default]")
    result.add_option("-r", unit=units.parsec, dest="rcl", type="float",
                      default=10 | units.parsec,
                      help="cluster half-mass radius [%default")
    result.add_option("-A", dest="algorithm", type="string",
                      default=BHTree, help="algorithm choice [%default]")
    #result.add_option("-a", dest="assignment", type="string",
    #                  default=None, help="Assignmentchoice [%default]")
    return result


if __name__ in '__main__':
    set_printing_strategy("custom",
                          preferred_units=[units.MSun, units.RSun, units.yr],
                          precision=4, prefix="", separator=" [",
                          suffix="]")
    o, arguments = new_option_parser().parse_args()
    nbody_integrator(**o.__dict__)

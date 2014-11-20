#!/usr/bin/env amuse
"""
SD_TLRH_s1603221.py

Computational Astrophysics / Stellar Evolution
Simulating Stars

Timo Halbesma, s1603221.
November 11, 2014. Version 1.0.
"""


import traceback

from amuse.lab import *
#from amuse.lab import EVtwin, SSE, SeBa
#from amuse.datamodel import Particle
#from amuse.community.mesa.interface import MESA
#from amuse.lab import Particles, new_salpeter_mass_distribution
#from amuse.lab import new_king_model, ph4
#from amuse.units import units, nbody_system


def evolve_a_single_star(M, z, model_time, algorithm="EVtwin()"):
    """
    @param M: mass of the star.
    @param z: metallicity of the star.
    @param stellar: specify the stellar evolution code.
        valid options are:
            MESA(), EVtwin() (both Henyey)
            SSE(), SeBa() (both parameterized)
    @return:

    """

    stellar = eval(algorithm)
    print "Algorithm used is {0}".format(algorithm)

    stellar.parameters.metallicity = z


    #stellar.commit_parameters()
    # print "\nstellar.parameters:\n{0}".format(stellar.parameters)
    star = Particles(mass=M)
    # print "\nStar:\n{0}".format(star)
    # print "\nstellar.particles:\n{0}".format(stellar.particles)
    stellar.particles.add_particles(star)

    channel_from_stellar = stellar.particles.new_channel_to(star)

    stellar.evolve_model(model_time)
    channel_from_stellar.copy()
    print "L(t=", star.age, ")=", star.luminosity, "R=", star.radius
    stellar.stop()


# 2D
def assignment_2d():
    N = 1024
    Rvir = 10 | units.parsec
    masses = new_salpeter_mass_distribution(N)
    converter = nbody_system.nbody_to_si(masses.sum(), Rvir)
    bodies = new_king_model(N, W0, convert_nbody=converter)
    bodies.mass = masses
    bodies.scale_to_standard(converter)

    # star the gravity solver
    gravity = ph4(converter)
    gravity.particles.add_particles(bodies)

    # star the stellar evolution solver
    stellar = MESA()
    stellar.particles.add_particles(bodies)

if __name__ == '__main__':
    # http://adsabs.harvard.edu/abs/2006CoAst.147...76A
    # metallicity = 0.0122

    # http://adsabs.harvard.edu/abs/2007ApJ...670..872C
    # metallicity = 0.0187 - 0.0239


    for metallicity in [0.02, 0.018, 0.0122, 0.0187, 0.0239]:
    # pass redirection=\"none" to receive full verbose output.
        for algorithm in ["SSE()", "SeBa()", "EVtwin()", "MESA()"]:
            try:
                evolve_a_single_star(1 | units.MSun, metallicity,
                                     4.6 | units.Gyr, algorithm)
            except Exception, e:  # filty generic, bad practice
                print Exception, e
                print traceback.format_exc()

            #raw_input("Press enter to continue")

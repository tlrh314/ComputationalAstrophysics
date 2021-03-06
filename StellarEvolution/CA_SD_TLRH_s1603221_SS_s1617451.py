#!/usr/bin/env amuse
"""
CA_SD_TLRH_s1603221_SS_s1617451.py

Computational Astrophysics / Stellar Evolution
Simulating Stars

Timo Halbesma, s1603221;
Shabaz Sultan, s1617451.

December 2, 2014. Version 1.0.
"""


import traceback
from time import time
import pickle

#from amuse.lab import *
from amuse import datamodel
from amuse.lab import EVtwin, SSE, SeBa
from amuse.community.mesa.interface import MESA

from amuse.lab import Particles

from amuse.lab import new_plummer_sphere, ph4
from amuse.units import units, nbody_system
from amuse.lab import set_printing_strategy
from amuse.community.bhtree.interface import BHTree

from amuse.rfi.core import is_mpd_running
from amuse.ic.salpeter import new_salpeter_mass_distribution


def evolve_a_single_star(m, z, model_time, algorithm=EVtwin):
    """
    Code for Assignment 2A.

    @param m: mass of the star.
    @param z: metallicity of the star.
    @param stellar: specify the stellar evolution code.
        valid options are:
            MESA, EVtwin (both Henyey)
            SSE, SeBa (both parameterized)
    @return: Dict with input parameters. mass/luminosity and runtime
    """

    t_start = time()

    stellar = algorithm(redirection="none")
    # stellar = algorithm()
    algorithm_name = str(algorithm).split('.')[-1][:-2]

    stellar.parameters.metallicity = z

    # stellar.commit_parameters()
    # print "\nstellar.parameters:\n{0}".format(stellar.parameters)
    star = Particles(mass=m)
    # print "\nStar:\n{0}".format(star)
    # print "\nstellar.particles:\n{0}".format(stellar.particles)
    stellar.particles.add_particles(star)

    channel_from_stellar = stellar.particles.new_channel_to(star)

    stellar.evolve_model(model_time)
    channel_from_stellar.copy()
    # print "L(t=", star.age, ")=", star.luminosity, "R=", star.radius, "\n"
    stellar.stop()

    runtime = (time() - t_start)

    set_printing_strategy("custom", preferred_units=[units.MSun, units.RSun,
                          units.Gyr, units.LSun], precision=4, prefix="",
                          separator=" [", suffix="]")
    # LaTeX table/tabular output <3
    print "{0} & {1} & {2} & {3} & {4} & {5} & {6} \\\\ "\
        .format(algorithm_name, z, model_time, star.luminosity,
                m, star.radius, runtime)

    data = {}
    data["algorithm_name"] = algorithm_name
    data["mass"] = m
    data["metallicity"] = z
    data["luminosity"] = star.luminosity
    data["radius"] = star.radius
    data["runtime"] = runtime

    #print data.items()

    return data


def create_cluster_isochrone(cluster, algorithm=SSE,
                             number_of_stars=1339,
                             end_time=1000 | units.Myr,
                             z=0.012
                             ):
    """
    Code for Assignment 2B.
    @param cluster: Name of cluter (used as dir to store data)
    @param stellar: specify the stellar evolution code.
        valid options are:
            MESA, EVtwin (both Henyey)
            SSE, SeBa (both parameterized)
    @param number_of_stars: how many stars to evolve
    @param end_time: how long the cluster is evolved
    @return: dumps the data to file
    """

    data = {}

    stellar = algorithm()
    algorithm_name = str(algorithm).split('.')[-1][:-2]

    print "Deriving a set of {0} random masses ".format(number_of_stars),
    print "following a Salpeter IMF between 0.1 and 125 Msun",
    print "(alpha = -2.35)."
    salpeter_masses = new_salpeter_mass_distribution(number_of_stars)

    print "Initializing the particles"
    stars = datamodel.Particles(number_of_stars)
    stars.mass = salpeter_masses

    total_mass = 0 | units.MSun
    for m in stars.mass: total_mass += m
    print "total mass =", total_mass

    data['initial_total_mass'] = total_mass

    stars.metallicity = z

    print "The evolution of {0} stars will be".format(number_of_stars),
    print "simulated until t = {0}".format(end_time),
    print ", using algorithm {0}".format(algorithm_name)
    stellar.commit_parameters()

    print "Stars to evolve:"
    print stars
    data['stars'] = stars

    stars = stellar.particles.add_particles(stars)
    stellar.commit_particles()

    print "Start evolving..."
    # stellar.evolve_model(end_time)

    times = [] | units.Myr
    luminosity_at_time = [] | units.LSun
    temperatures_at_time = [] | units.K
    mass_at_time = [] | units.MSun

    current_time = 0 | units.Myr
    while current_time < end_time:
        name_of_the_figure = "isochrone_{0}_{1}_{2}.png".\
            format(algorithm_name, number_of_stars,
                   int(current_time.value_in(units.Myr)))

        stellar.evolve_model(current_time)
        temperatures = stars.temperature
        luminosities = stars.luminosity

        times.append(current_time)
        luminosity_at_time.append(luminosities)
        temperatures_at_time.append(temperatures)
        mass_at_time.append(stars.total_mass())


        print current_time
        plot_HR_diagram(temperatures, luminosities,
                        cluster+"/",
                        name_of_the_figure, current_time)

        current_time += 10 | units.Myr

    print "Evolved model successfully."
    stellar.stop()

    print "All done!"

    data['times'] = times
    data['luminosity_at_time'] = luminosity_at_time
    data['temperatures_at_time'] = temperatures_at_time
    data['mass_at_time'] = mass_at_time

    print data
    pickle.dump(data, open(cluster+"/assignment2b.dat", "wb"))


def plot_HR_diagram(temperatures, luminosities, directory,
                    name_of_the_figure, end_time):
    try:
        import matplotlib
        # This removes the need for ssh -X to enable plotting.
        matplotlib.use("Agg", warn=False)
        from matplotlib import pyplot
    except ImportError:
        print "Unable to produce plot: couldn't find matplotlib."
    else:
        print "Plotting the data..."
        number_of_stars = len(temperatures)
        pyplot.figure(figsize=(7, 8))
        pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
        pyplot.xlabel(r'T$_{\rm eff}$ (K)')
        pyplot.ylabel(r'Luminosity (L$_\odot$)')
        pyplot.loglog(temperatures.value_in(units.K),
                      luminosities.value_in(units.LSun), "ro")

        xmin, xmax = 20000.0, 2500.0
        ymin, ymax = 1.e-4, 1.e4
        pyplot.text(xmin*.75, ymax*0.1, str(number_of_stars) +
                    " stars\nt="+str(end_time))
        pyplot.axis([xmin, xmax, ymin, ymax])
        pyplot.savefig(directory+name_of_the_figure)
        pyplot.close()


def assignment_2a():
    # http://adsabs.harvard.edu/abs/2006CoAst.147...76A
    # metallicity = 0.0122

    # http://adsabs.harvard.edu/abs/2007ApJ...670..872C
    # metallicity = 0.0187 - 0.0239

    datapoints = []
    for algorithm in [MESA, EVtwin, SSE, SeBa]:
    # pass redirection=\"none" to receive full verbose output.
        for metallicity in [0.02]: # [0.02, 0.018, 0.0122, 0.0187, 0.0239]:
            try:
                datapoints.append(evolve_a_single_star(1 | units.MSun,
                                  metallicity, 4.6 | units.Gyr, algorithm))
                datapoints.append(evolve_a_single_star(10 | units.MSun,
                                  metallicity, 30 | units.Myr, algorithm))
                # datapoints.append(evolve_a_single_star(100 | units.MSun,
                #                   metallicity, 4.6 | units.Myr, algorithm))
            except Exception, e:  # filty generic, bad practice
                print Exception, e
                print traceback.format_exc()

            # raw_input("Press enter to continue")
    pickle.dump(datapoints, open("assignment2a.dat", "wb"))


def assignment_2b():
    # Metallicity according to WEBDA page for open cluster Hyades
    create_cluster_isochrone("Hyades", SSE, 1339,
                             1000 | units.Myr, 0.17)

    # create_cluster_isochrone("Pleiades", SSE, 887,
    #                          1000 | units.Myr, 0.17)


def stellar_remnant_state(star):
    return 10 <= star.stellar_type.value_in(units.stellar_type) and \
        star.stellar_type.value_in(units.stellar_type) < 16

def stellar_main_sequence_state(star):
    return star.stellar_type.value_in(units.stellar_type) in (0,1,7)

def stellar_giant_state(star):
    return star.stellar_type.value_in(units.stellar_type) in (2,3,4,5,6,8,9)

def assignment_2d():
    current_cluster_mass = 400 | units.MSun
    initial_mass_fraction = 0.84
    desired_initial_mass = current_cluster_mass / initial_mass_fraction

    masses = new_salpeter_mass_distribution(100000)
    mean_salpeter_mass = masses.mean()
    print "mean salpeter mass", mean_salpeter_mass
    N = int(desired_initial_mass / mean_salpeter_mass)
    print "N", N


    Rvir = 10 | units.lightyear
    z = 0.17
    masses = new_salpeter_mass_distribution(N)
    converter = nbody_system.nbody_to_si(masses.sum(), Rvir)
    G_SI = converter.to_si(nbody_system.G)
    bodies = new_plummer_sphere(N, convert_nbody=converter)
    bodies.mass = masses
    bodies.metalicity = z

    # start the gravity solver
    gravity = BHTree(converter)
    gravity.initialize_code()
    gravity.parameters.timestep = 0.1 | units.Myr

    # start the stellar evolution solver
    stellar = SSE()
    stars = stellar.particles.add_particles(bodies)
    from_stellar_evolution_to_model \
        = stellar.particles.new_channel_to(bodies)
    from_stellar_evolution_to_model.copy_attributes(["mass"])

    bodies.scale_to_standard(converter)
    gravity.particles.add_particles(bodies)

    from_model_to_gravity = bodies.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(bodies)
    gravity.commit_particles()

    end_time = 1000 | units.Myr
    current_time = 0 | units.Myr
    cluster = "Hyades"
    bound_stars_counts = []
    main_sequence_stars_counts = []
    giant_stars_counts = []
    remnant_stars_counts = []
    max_radii = [] | units.parsec
    virial_radii = [] | units.parsec
    times = [] | units.Myr
    while current_time < end_time:
        name_of_the_figure = "isochrone_with_grav_"+str(int(current_time.value_in(units.Myr)))+".png"

        gravity.evolve_model(current_time)
        stellar.evolve_model(current_time)


        from_gravity_to_model.copy()
        from_stellar_evolution_to_model.copy_attributes(["mass", "radius"])
        from_model_to_gravity.copy_attributes(["mass"])


        remnant_count = 0
        main_sequence_count = 0
        giant_count = 0
        for star in stars:
            if stellar_remnant_state(star):
                remnant_count += 1
            if stellar_giant_state(star):
                giant_count += 1
            if stellar_main_sequence_state(star):
                main_sequence_count += 1
        max_radius = bodies.total_radius()
        virial_radius = bodies.virial_radius()
        bound_star_count = len(bodies.bound_subset(unit_converter=converter, G=G_SI))
        print "bound stars:", bound_star_count
        print "main sequence stars:", main_sequence_count
        print "giant stars:", giant_count
        print "remnant stars:", remnant_count
        print "cluster radius(max from centre):", max_radius
        print "virial radius:", virial_radius
        print current_time

        times.append(current_time)
        remnant_stars_counts.append(remnant_count)
        giant_stars_counts.append(giant_count)
        main_sequence_stars_counts.append(main_sequence_count)
        max_radii.append(max_radius)
        virial_radii.append(virial_radius)
        bound_stars_counts.append(bound_star_count)


        temperatures = stars.temperature
        luminosities = stars.luminosity

        plot_HR_diagram(temperatures, luminosities,
                        cluster+"/",
                        name_of_the_figure, current_time)
        current_time += 10 | units.Myr
    data = {}
    data["bound_stars_at_time"] = bound_stars_counts
    data["remnant_stars_at_time"] = remnant_stars_counts
    data["giant_stars_at_time"] = giant_stars_counts
    data["main_sequence_stars_at_time"] = main_sequence_stars_counts
    data["max_radius_at_time"] = max_radii
    data["virial_radii"] = virial_radii
    data["times"] = times
    pickle.dump(data, open(cluster+"/assignment2d.dat", "wb"))


if __name__ == '__main__':
    assert is_mpd_running()
    #assignment_2a()
    #assignment_2b()
    assignment_2d()

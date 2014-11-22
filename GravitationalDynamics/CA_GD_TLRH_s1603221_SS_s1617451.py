#!/usr/bin/env amuse
"""
GD_TLRH_s1603221.py

Computational Astrophysics / Gravitational Dynamics
Simulating Ultra Compact Dwarf Galaxies with AMUSE

Timo Halbesma, s1603221.
Shebaz Sultan, s1617451.

November 25, 2014. Version 1.0.
"""

from time import time

import numpy
from matplotlib import pyplot
from amuse.units import units
from amuse.community.hermite0.interface import Hermite
from amuse.community.huayno.interface import Huayno
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE

from solve_nbody import nbody_integrator


def assignment_1a(algorithm=BHTree):
    options = dict()
    stars = numpy.array([1, 2, 4, 8, 16, 32, 64, 128, 256,
                         512, 1024], dtype=numpy.int16)
    t_end = numpy.array([1, 2, 4, 8, 16, 32], dtype=numpy.int8)

    options['mcl'] = 10**7 | units.MSun
    options['rcl'] = 10 | units.parsec
    options['n_steps'] = 100
    options['algorithm'] = algorithm

    runtime_of_N_and_t = numpy.zeros((len(stars), len(t_end)),
                                     dtype=numpy.float64)
    dE_of_N_and_t = numpy.zeros((len(stars), len(t_end)),
                                dtype=numpy.float64)

    for i, N in enumerate(stars):
        for j, t in enumerate(t_end):
            options['Ncl'] = N
            options['t_end'] = t | units.Myr
            t_start = time()
            dE_of_N_and_t[i][j] = nbody_integrator(**options)
            runtime_of_N_and_t[i][j] = (time() - t_start)
            print str(N) + ', ' + str(t) + ', ' + str(time() - t_start)\
                + ', ' + str(dE_of_N_and_t[i][j])
    return stars, t_end, dE_of_N_and_t, runtime_of_N_and_t,\
        options['algorithm']


def plot_1a(stars, t_end, to_plot, choice, algorithm):
    colors = {0: 'r', 1: 'y', 2: 'g', 3: 'b', 4: 'c', 5: 'm'}
    labelled = False
    fig, ax = pyplot.subplots()
    algorithm_name = str(algorithm).split('.')[-1][:-2]

    for i, N in enumerate(stars):
        for j, integration_time in enumerate(t_end):
            print "N = {0}, t = {1}, dE = {2}, to_plot = {3}"\
                  .format(N, integration_time, to_plot[i][j], dE[i][j])
            if not labelled:
                ax.scatter(N, to_plot[i][j],
                           color=colors.get(j, 'k'),
                           label="t={0} Myr".format(integration_time))
            else:
                ax.scatter(N, to_plot[i][j],
                           color=colors.get(j, 'k'))
        labelled = True

    if choice == "runtime":
        ax.set_title("Wall-clock runtime as a function of N and t_end " +
                     "using algorithm '{0}'".format(algorithm_name))
        ax.set_ylabel("Wall-clock time (s)")
    elif choice == "dE":
        ax.set_title("Relative energy error as a function of N and t_end " +
                     "using algorithm '{0}'".format(algorithm_name))
        ax.set_ylabel("relative energy error dE")

    ax.set_xlabel("N")
    ax.set_xscale('log', basex=2)
    ax.legend(loc=2)  # upper left
    pyplot.savefig("plots/CA_GD_TLRH_s1603221_SS_s1617451_{0}_{1}"
                   .format(algorithm_name, choice))


def radius_dependency(algorithm, radii, N, t):
    options = dict()

    runtime_of_r = numpy.zeros(len(radii), dtype=numpy.float64)

    options['Ncl'] = N
    options['mcl'] = 10**7 | units.MSun
    options['t_end'] = t | units.Myr
    options['n_steps'] = 100
    options['algorithm'] = algorithm

    for i, radius in enumerate(radii):
        options['rcl'] = radius | units.parsec
        t_start = time()
        nbody_integrator(**options)
        runtime_of_r[i] = (time() - t_start)
        print "Runtime: {0}".format(runtime_of_r[i])

    return runtime_of_r


def plot_radius_dependency(radii, N, t, runtimes):
    colors = {'BHTree': 'r', 'Hermite': 'y', 'Huayno': 'g', 'PhiGRAPE': 'b',
              'Mercury': 'c', 'mmc': 'm'}

    fig, ax = pyplot.subplots()

    for algorithm, times in runtimes.items():
        ax.scatter(radii, times,
                   color=colors.get(algorithm, 'k'),
                   label=algorithm)

    ax.set_title("Radius Dependency for N={0}, t_end={1}"
                 .format(N, t))
    ax.set_xlabel("Cluster Half-Mass Radius rcl (parsec)")
    ax.set_ylabel("Wall-clock time (s)")

    ax.set_xscale('log', basex=2)
    ax.legend(loc=2)  # upper left
    pyplot.savefig("plots/CA_GD_TLRH_s1603221_SS_s1617451_rcldep_N{0}_t-end{1}"
                   .format(N, t))


if __name__ in '__main__':
    # Assignment 1A
    # stars, t_end, dE, runtime, algorithm = assignment_1a(Hermite)
    # plot_1a(stars, t_end, runtime, "runtime", algorithm)
    # plot_1a(stars, t_end, dE, "dE", algorithm)
    # pyplot.show()  # Only works in AMUSE 8.1 (binary)

    # Assignment 1B
    # To implement...

    # Assignment 1C
    runtime_of_radii = dict()
    radii = numpy.array([1, 2, 4, 8, 16, 32, 64, 128, 256,
                         512, 1024], dtype=numpy.int16)

    for N in [8, 16, 32]:
        for t in [8, 16, 32]:
            for integrator in [BHTree, Hermite, Huayno]:
                # stars, t_end, dE, runtime, algorithm = assignment_1a(integrator)
                # plot_1a(stars, t_end, runtime, "runtime", algorithm)
                # plot_1a(stars, t_end, dE, "dE", algorithm)
                # pyplot.show()  # Only works in AMUSE 8.1 (binary)

                algorithm_name = str(integrator).split('.')[-1][:-2]
                runtime_of_radii[algorithm_name] = radius_dependency(integrator, radii, N, t)
                # raw_input("Press enter to continue")

            plot_radius_dependency(radii, N, t, runtime_of_radii)

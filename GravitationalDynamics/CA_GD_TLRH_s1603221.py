#!/usr/bin/env amuse
"""
computationalastrophysics_gravitationaldynamics_timohalbesma_s1603221.py

Computational Astrophysics / Gravitational Dynamics
Simulating Ultra Compact Dwarf Galaxies with AMUSE / Excersize 1A.

Timo Halbesma, s1603221.
October 14, 2014. Version 1.0.
"""
from time import time

from amuse.units import units

from solve_nbody import nbody_integrator


if __name__ in '__main__':
    options = dict()
    stars = [1, 2, 4, 8, 16, 32, 65, 128, 256, 512, 1024]
    dE_of_N_and_t = list()
    runtime_of_N_and_t = list()
    t_end = [1, 2, 4, 8, 16, 32]

    options['mcl'] = 10**7 | units.MSun
    options['rcl'] = 10 | units.parsec
    options['n_steps'] = 100
    options['algorithm'] = "BHTree"

    for N in stars:
        dE_of_t = list()
        runtime_of_t = list()
        for t in t_end:
            options['Ncl'] = N
            options['t_end'] = t | units.Myr
            # print "N = {0}, t = {1}, options = {2}".format(N, t, options)
            t_start = time()
            dE_of_t.append(nbody_integrator(**options))
            runtime_of_t.append((time() - t_start) | units.s)
        dE_of_N_and_t.append(dE_of_t)
        runtime_of_N_and_t.append(runtime_of_t)
        break

    print dE_of_N_and_t
    print runtime_of_N_and_t

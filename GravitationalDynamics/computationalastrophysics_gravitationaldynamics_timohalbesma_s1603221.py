#!/usr/bin/env amuse
"""
computationalastrophysics_gravitationaldynamics_timohalbesma_s1603221.py

Computational Astrophysics / Gravitational Dynamics
Simulating Ultra Compact Dwarf Galaxies with AMUSE / Excersize 1A.

Timo Halbesma, s1603221.
October 14, 2014. Version 1.0.
"""
import cProfile

import numpy as np

from solve_nbody import nbody_integrator, new_option_parser


if __name__ in '__main__':
    print np.array(1, 2, 4, 8, 16, 32, 65, 128, 256, 512, 1024)
    print np.array(1, 2, 4, 8, 16, 32)
    # options, arguments = new_option_parser().parse_args()

    # cProfile.run('''nbody_integrator(**options.__dict__)''')

import numpy
from matplotlib import pyplot

from amuse.community.bhtree.interface import BHTree
from amuse.units import units

from solve_nbody import nbody_integrator

# Assignment 1B
if __name__ in '__main__':
    options = {}
    options['mcl'] = 10**7 | units.MSun
    options['rcl'] = 10 | units.parsec
    options['n_steps'] = 100
    options['algorithm'] = BHTree
    N = 1024
    options['Ncl'] = N
    t_end = 2 | units.Myr
    options['t_end'] = t_end
    error_for_timestep = []
    trials = 10
    for i in xrange(11):
        timestep = 1.0/2**i
        options['timestep'] = timestep
        error = 0
        for j in xrange(trials):
            error += abs(nbody_integrator(**options))
        print timestep, error
        error_for_timestep.append((timestep, error/trials))

    pyplot.scatter([i[0] for i in error_for_timestep],
            [i[1] for i in error_for_timestep])
    pyplot.xlabel("timestep (fraction of dt)")
    pyplot.ylabel("magnitude of relative error")
    pyplot.xscale("log", basex=2)
    pyplot.xlim(0,1)
    pyplot.xticks([i[0] for i in error_for_timestep][::-1])
    pyplot.title(
            "Relative Energy Error as function of timestep for Barnes-Hut, N="
            +str(N)+",\n t_end="+str(t_end)+", "+str(trials)+
            " trials per datapoint")
    pyplot.savefig("plots/error_for_timestep")

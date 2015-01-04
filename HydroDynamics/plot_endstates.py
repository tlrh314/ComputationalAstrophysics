import os

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

from amuse import plot as aplot
from amuse.io import read_set_from_file
from amuse.units import units

if __name__ == "__main__":
    # for i in xrange(5, 21, 5):
    # for i in xrange(20, 301, 20):
    for i in [20]:
        impact_parameter = i | units.parsec
        # for j in xrange(5, 26, 5):
        # for j in xrange(25, 201, 25):
        for j in [15]:
            v_infinity = (j/100.0) | units.kms

            data_file = "simulation_data/v_"+str(j)+"_d_"+str(i)+".dat"
            out_file = "simulation_end_states/v_"+str(j)+"_d_"+str(i)+".png"

            # if os.path.isfile(out_file):
            #    print "file " + out_file + " already exists, skip plotting"
            #    continue
            print "ouputing " + out_file
            cluster, gascloud = read_set_from_file(data_file, 'amuse',
                                                   names=("cluster", "gas"))
            lim = abs(next(cluster.history).position).max()

            fig = pyplot.figure(figsize=[10, 10])

            artists = []
            cl = None
            gas = None
            for c, g in zip(cluster.history, gascloud.history):
                cl = c
                gas = g
            points = aplot.scatter(cl.x, cl.y, c='black',
                                   label='Stellar Cluster')
            gas_points = aplot.scatter(gas.x, gas.y, s=100, c='b',
                                       edgecolors='none', alpha=0.3,
                                       rasterized=True,
                                       label='Giant Molecular Cloud')

            time_label = r'$t_{\rm end}$ ' + "= {:.2f} Myr".format(cl.get_timestamp().
                                                 value_in(units.Myr))
            text = aplot.text(-4./5.*lim, 4./5.*lim, time_label)
            aplot.xlabel("X")
            aplot.ylabel("Y")
            pyplot.axis('equal')
            aplot.xlim(-lim, lim)
            aplot.ylim(-lim, lim)
            pyplot.legend()
            pyplot.title(r'Cluster orbitting GMC with $d$={0}, $v_\infty$={1}'
                         .format(impact_parameter, v_infinity) +
                         r', at $t_{\rm end}$')
            pyplot.savefig(out_file)
            pyplot.close()

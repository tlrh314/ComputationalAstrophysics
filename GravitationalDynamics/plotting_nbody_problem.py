"""
    Visualization for simple N-body integration.
    Reads particle set from file (nbody.hdf5) and prints
        subsequent frames.
"""
from matplotlib import pyplot
from amuse.lab import *
from amuse.plot import scatter, xlabel, ylabel


def plot_cluster(filename="nbody.hdf5"):
    """ Plot file nbody.hdf5 """
    pyplot.ion()
    stars = read_set_from_file(filename, format='hdf5')
    lim = 10 * stars.center_of_mass().length().value_in(stars.x.unit)
    m = 1 + 3.0 * stars.mass / min(stars.mass)
    for si in stars.history:
        time = si.get_timestamp()
        pyplot.title("Cluster at t="+str(time))
        print "time =", time
        scatter(si.x, si.y, s=m)
        xlabel("X")
        ylabel("Y")
        pyplot.xlim(-lim, lim)
        pyplot.ylim(-lim, lim)
        pyplot.draw()
        pyplot.cla()
    pyplot.show()


def new_option_parser():
    """ Set options """
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default="nbody.hdf5",
                      help="output filename [nbody.hdf5]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments = new_option_parser().parse_args()
    plot_cluster(**o.__dict__)

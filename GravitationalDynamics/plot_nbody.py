"""
    Visualization for simple N-body integration.
    Reads particle set from file (nbody.hdf5) and prints
        subsequent frames.
"""
from matplotlib import pyplot
from amuse.lab import read_set_from_file  # take a hike with your *
from amuse import plot as aplot # scatter, xlabel, ylabel


def plot_cluster(filename="nbody.hdf5"):
    """ Plot file nbody.hdf5 """
    pyplot.ion()  # Turn interactive mode on.
    stars = read_set_from_file(filename, format='hdf5')
    lim = 1000*stars.center_of_mass().length().value_in(stars.x.unit)
    m = 1 + 3.0*stars.mass/min(stars.mass)
    print "m =", m

    for si in stars.history:
        time = si.get_timestamp()
        pyplot.title("Cluster at t="+str(time))
        print "time =", time
        aplot.scatter(si.x, si.y, s=m)  # s size (in point^2)
        aplot.xlabel("X")
        aplot.ylabel("Y")
        pyplot.xlim(-lim, lim)
        pyplot.ylim(-lim, lim)
        pyplot.draw()  # Redraw the current figure (interactive mode).
        pyplot.cla()  # Clear the current axes.


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

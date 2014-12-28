import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot, animation
from amuse import plot as aplot

from amuse.support.console import set_printing_strategy
from amuse.io import read_set_from_file
from amuse.units import units, optparse

def make_movie(input_filename, output_filename):
    cluster, gas = read_set_from_file(input_filename, 'amuse', 
            names=("cluster", "gas"))
    lim = abs(next(cluster.history).position).max()
    
    fig = pyplot.figure(figsize=[10,10])

    artists = []
    for cl, gas in zip(cluster.history, gas.history):
        print "Creating frame at time", cl.get_timestamp()
        points = aplot.scatter(cl.x, cl.y, c='black')
        gas_points = aplot.scatter(gas.x, gas.y, s=100, c='b', 
                edgecolors='none', alpha=0.3, rasterized=True)

        time_label = "t = {:.2f} Myr".format(
                cl.get_timestamp().value_in(units.Myr))
        text = aplot.text(-4./5.*lim, 4./5.*lim, time_label)
        aplot.xlabel("X")
        aplot.ylabel("Y")
        pyplot.axis('equal')
        aplot.xlim(-lim, lim)
        aplot.ylim(-lim, lim)
        pyplot.savefig("plots/test."+str(cl.get_timestamp().value_in(units.Myr))+".png")
        pyplot.close()
                

        artists.append((points, gas_points, text))

    #aplot.xlabel("X")
    #aplot.ylabel("Y")
    #pyplot.axis('equal')
    #aplot.xlim(-lim, lim)
    #aplot.ylim(-lim, lim)

    #print "Writing movie to file"
    #movie = animation.ArtistAnimation(fig, artists, interval=50,
    #        repeat_delay=3000, blit=True)
    #writer = animation.writers['ffmpeg'](fps=15, metadata=dict(artist='Me'),
    #    bitrate=1800)
    #movie.save(output_filename, writer=writer, dpi=300)

def parse_options():
    parser= optparse.OptionParser()
    parser.add_option("-i", dest="input_filename", 
            type="string", default="test.amuse", 
            help="The file in which the snapshots are saved [%default]")
    parser.add_option("-o", dest="output_filename",
            type="string", default="test.mp4",
            help="The name of the movie [%default]")

    options, arguments = parser.parse_args()
    return options.__dict__

if __name__ == "__main__":
    set_printing_strategy("custom", preferred_units = [units.parsec, units.Myr])
    options = parse_options()
    make_movie(**options)

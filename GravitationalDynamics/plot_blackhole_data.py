from matplotlib import pyplot
import pickle 


from amuse.plot import scatter, xlabel, ylabel

if __name__ in '__main__':
    datapoints = pickle.load(open("bound_mass.dat", "rb"))
    filename = "plots/bound_mass_at_6000_pc"

    pyplot.title("bound mass fraction of cluster orbiting at 6kpc")
    for datapoint in datapoints:
        scatter(datapoint["radius"], datapoint["hop_mass_fraction"])
    xlabel("cluster radius")
    ylabel("bound fraction of cluster mass")
    pyplot.savefig(filename)

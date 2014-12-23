import pickle

from amuse.units import units

from matplotlib import pyplot
from amuse.plot import xlabel, ylabel, plot

data = pickle.load(open("assignment2d.dat"))


pyplot.figure(figsize = (8, 6))
plot(data["times"], data["giant_stars_at_time"], "g", label="Giant Stars")
plot(data["times"], data["main_sequence_stars_at_time"], "r", label="Main Sequence Stars")
plot(data["times"], data["remnant_stars_at_time"], "b", label="Remnant Stars")
pyplot.legend()
pyplot.title("stars in different stellar phases over time")
pyplot.savefig("stellar_phases_counts")
pyplot.close()

pyplot.figure(figsize = (8, 6))
plot(data["times"], data["giant_stars_at_time"], "g", label="Giant Stars")
plot(data["times"], data["remnant_stars_at_time"], "b", label="Remnant Stars")
pyplot.legend()
pyplot.title("stars in different stellar phases over time")
pyplot.savefig("stellar_phases_counts_without_main_sequence")
pyplot.close()

pyplot.figure(figsize = (8, 6))
plot(data["times"], data["bound_stars_at_time"])
pyplot.title("number of bound stars")
pyplot.savefig("bound_stars")
pyplot.close()

pyplot.figure(figsize = (8, 6))
plot(data["times"], data["max_radius_at_time"])
pyplot.title("radius of cluster (max distance of a star from centre)")
pyplot.savefig("cluster_max_radius")
pyplot.close()

pyplot.figure(figsize = (8, 6))
plot(data["times"], data["virial_radii"])
pyplot.title("virial radius of cluster")
pyplot.savefig("cluster_virial_radius")
pyplot.close()

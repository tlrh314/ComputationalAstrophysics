import pickle

from amuse.units import units

from matplotlib import pyplot
from amuse.plot import xlabel, ylabel, plot


data = pickle.load(open("assignment2b.dat"))

masses = data['mass_at_time']
initial_mass = data['initial_total_mass']
times = data['times']
mass_fractions = [m/initial_mass for m in masses]

pyplot.figure(figsize = (8, 6))
plot(times, masses)
pyplot.title("mass in cluster")
pyplot.savefig("mass_over_time")
pyplot.close()

pyplot.figure(figsize = (8, 6))
plot(times, mass_fractions)
pyplot.title("fraction of initial mass in cluster")
pyplot.ylabel("mass fraction")
pyplot.savefig("massfraction_over_time")
pyplot.close()

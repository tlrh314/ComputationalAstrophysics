from amuse.units import units

from matplotlib import pyplot
from amuse.plot import xlabel, ylabel


f = open("fitting_data.dat")
times = []
fitness = []

for line in f.readlines():
    words = line.split()
    fitness.append(float(words[0]))
    times.append(float(words[1]))

pyplot.figure(figsize = (8, 6))
pyplot.plot(times, fitness)
pyplot.title("squared difference with simulation data over time for Hyades cluster")
pyplot.xlabel("time in Myr")
pyplot.ylabel("squared difference with sim data")
pyplot.savefig("fitness_over_time")
pyplot.close()

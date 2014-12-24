import pickle
import math
from amuse.units import units

from matplotlib import pyplot
from amuse.plot import xlabel, ylabel

import sys

data = pickle.load(open("assignment2b.dat"))

luminosity_at_time = data['luminosity_at_time']
total_mass = data['total_mass']
temperatures_at_time = data['temperatures_at_time']
stars = data['stars']
times = data['times']

average_luminosity = 0
average_temperature = 14
simulated_star_count = 0
for timestep in luminosity_at_time:
    for luminosity in timestep:
        average_luminosity += luminosity.number
        simulated_star_count += 1
for timestep in temperatures_at_time:
    for temperature in timestep:
        average_temperature += temperature.number

print "based on sim data:"
average_luminosity /= simulated_star_count
average_temperature /= simulated_star_count
print "average luminosity", average_luminosity
print "average temperature", average_temperature

luminosity_std_dev = 0
for timestep in luminosity_at_time:
    for luminosity in timestep:
        luminosity_std_dev += (luminosity.number - average_luminosity)**2
luminosity_std_dev = math.sqrt(luminosity_std_dev/simulated_star_count)

temperature_std_dev = 0
for timestep in temperatures_at_time:
    for temperature in timestep:
        temperature_std_dev += (temperature.number - average_temperature)**2
temperature_std_dev = math.sqrt(temperature_std_dev/simulated_star_count)
print "luminosity std_dev:", luminosity_std_dev
print "temperature std_dev:", temperature_std_dev

f = open("ubv_hyades.peo")
f.readline()
f.readline()

#Perryman et al., 1998
hyades_distance = 46.34 | units.parsec

observational_data = []

observed_temps = []
observed_lums = []

for line in f.readlines():
    words = line.split()
    if len(words) != 6:
        words =  line.split('\t')
    datum = {}
    datum['star_number'] = int(words[0])
    datum['data_source'] = int(words[1])
    if words[2]:
        datum['V'] = float(words[2])
        Vmag = datum['V']
        solar_luminosity = (10**((Vmag + 26.7)/-2.5))/(((1|units.AU)/hyades_distance)**2)
        #solar_luminosity = hyades_distance**2 *  10**((Vmag+2.72)/-2.5)
        print solar_luminosity
        datum['solar_luminosity'] = solar_luminosity | units.LSun
        observed_lums.append(solar_luminosity)
    else:
        continue
    datum['B-V'] = float(words[3])
    bv = datum['B-V']
    # Ballesteros' formula
    # used http://en.wikipedia.org/wiki/Color_index as ref
    temp = 4600*(1.0/(0.92 * bv + 1.7) + 1.0/(0.92 * bv + 0.62))
    observed_temps.append(temp)
    datum['temperature'] = temp | units.K
    if words[4]:
        datum['U-B'] = float(words[4])
    datum['N'] = int(words[5])
    observational_data.append(datum)


#for datum in observational_data[:100]:
    #print datum['solar_luminosity']

def normalized_temperature_diff(temp_a, temp_b, average_temp, std_dev_temp):
    normalized_a = (temp_a.number - average_temp) / std_dev_temp
    normalized_b = (temp_b.number - average_temp) / std_dev_temp
    return normalized_b-normalized_a

def normalized_luminosity_diff(lum_a, lum_b, average_luminosity, std_dev_lum):
    normalized_a = (lum_a.number - average_luminosity) / std_dev_lum
    normalized_b = (lum_b.number - average_luminosity) / std_dev_lum
    return normalized_b-normalized_a

def plot_figure(filename, temperatures, luminosities):
        pyplot.figure(figsize = (7, 8))
        pyplot.loglog(temperatures.value_in(units.K), luminosities.value_in(units.LSun), "ro")
        #xmin, xmax = 20000.0, 2500.0
        #ymin, ymax = 1.e-4, 1.e4
        #pyplot.text(xmin*.75,ymax*0.1,str(number_of_stars)+" stars\nt="+str(end_time))
        #pyplot.axis([xmin, xmax, ymin, ymax])
        pyplot.savefig(filename)

for i, time in enumerate(times):
    squares_sum = 0
    for observed_star in observational_data:
        #find closest star in simulated dataset
        lowest_temp_diff = float('inf')
        lowest_luminosity_diff = float('inf')
        lowest_diff = float('inf')


        for simulated_star in zip(temperatures_at_time[i][:], luminosity_at_time[i][:]):
            temperature_diff = normalized_temperature_diff(observed_star['temperature'], simulated_star[0], average_temperature, temperature_std_dev)**2
            luminosity_diff = normalized_luminosity_diff(observed_star['solar_luminosity'], simulated_star[1], average_luminosity, luminosity_std_dev)**2
            difference = luminosity_diff + temperature_diff
            if difference < lowest_diff:
                lowest_diff = difference

        squares_sum += lowest_diff

    pyplot.figure(figsize = (7, 8))
    xmin, xmax = 20000.0, 2500.0
    ymin, ymax = 1.e-4, 1.e4
    pyplot.axis([xmin, xmax, ymin, ymax])
    pyplot.loglog(temperatures_at_time[i][:].value_in(units.K), luminosity_at_time[i][:].value_in(units.LSun), "g.", label="Simulated")
    observed_temperatures = [j["temperature"].number for j in observational_data]
    a = observed_temperatures | units.K
    observed_luminosities = [j["solar_luminosity"].number for j in observational_data]
    b = observed_luminosities | units.LSun
    pyplot.loglog(observed_temps, observed_lums, 'r.', label="Observed")
    pyplot.xlabel(r'$\log \, T_{\rm eff} (K)$')
    pyplot.ylabel(r'$\log \, L (L_{\odot})$')
    pyplot.title("HR diagram at time " + str(time) + r' with $\chi^2$=' + str(squares_sum))
    pyplot.legend()
    pyplot.savefig("test"+str(int(time.number)))
    pyplot.close()
    print squares_sum, time

    print i, time
    import sys; sys.exit()

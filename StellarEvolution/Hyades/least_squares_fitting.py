import pickle
import math
from amuse.units import units

data = pickle.load(open("assignment2b.dat"))

luminosity_at_time = data['luminosity_at_time']
total_mass = data['total_mass']
temperatures_at_time = data['temperatures_at_time']
print temperatures_at_time[0][:10]
print luminosity_at_time[0][:10]
stars = data['stars']
times = data['times']

f = open("ubv_hyades.peo")
f.readline()
f.readline()

hyades_distance = 146.7 | units.lightyear

observational_data = []

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
        solar_luminosity = (10**(Vmag - 26.7))/(((1|units.AU)/hyades_distance)**2)
        datum['solar_luminosity'] = solar_luminosity | units.LSun
    else:
        continue
    datum['B-V'] = float(words[3])
    bv = datum['B-V']
    # Ballesteros' formula
    # used http://en.wikipedia.org/wiki/Color_index as ref
    temp = 4600*(1.0/(0.92 * bv + 1.7) + 1.0/(0.92 * bv + 0.62))
    datum['temperature'] = temp | units.K
    if words[4]:
        datum['U-B'] = float(words[4])
    datum['N'] = int(words[5])
    observational_data.append(datum)

#for datum in observational_data[:100]:
    #print datum['solar_luminosity']


for i, time in enumerate(times):
    squares_sum = 0
    for observed_star in observational_data:
        #find closest star in simulated dataset
        lowest_temp_diff = float('inf') 
        lowest_luminosity_diff = float('inf')
        lowest_diff = float('inf')
        for simulated_star in zip(temperatures_at_time[i][::5], luminosity_at_time[i][::5]):
            temperature_diff = ((observed_star['temperature'] - simulated_star[0])**2).number
            if temperature_diff < lowest_temp_diff:
                lowest_temp_diff = temperature_diff
            luminosity_diff = ((observed_star['solar_luminosity'] - simulated_star[1])**2).number * 10**8
            if luminosity_diff < lowest_luminosity_diff:
                lowest_luminosity_diff = luminosity_diff
            difference = luminosity_diff + temperature_diff
            if difference < lowest_diff:
                lowest_diff = difference

            #print simulated_star
            #print observed_star['temperature']
            #print observed_star['solar_luminosity']
        #print "diff", lowest_diff
        #print "temp diff", lowest_temp_diff
        #print "luminosity diff", lowest_luminosity_diff
        squares_sum += lowest_diff
    print time
    print squares_sum
    #print i, time

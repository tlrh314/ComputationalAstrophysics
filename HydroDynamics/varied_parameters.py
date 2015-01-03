import os

from amuse.support.console import set_printing_strategy
from amuse.units import units

from CA_HD import evolve_cluster_and_cloud

if __name__ == "__main__":
    set_printing_strategy("custom",
            preferred_units = [units.MSun, units.parsec, units.Myr, units.kms])

    cluster_stars = 128
    cluster_mass = 10**6|units.MSun 
    cluster_radius = 10|units.parsec
    gas_particles = 1024
    gas_mass = 10**7|units.MSun
    gas_radius = 100|units.parsec
    impact_parameter = 150|units.parsec
    v_infinity = 1|units.kms
    timestep = 0.2|units.Myr
    filename = "test.amuse"

    
    for i in xrange(5, 21, 5):
        impact_parameter = i|units.parsec
        for j in xrange(5, 26, 5):
            v_infinity = (j/100.0)|units.kms

            filename = "simulation_data/v_"+str(j)+"_d_"+str(i)+".dat"

            exists = os.path.isfile(filename)

            if not exists:
                print "creating " + filename
                evolve_cluster_and_cloud(cluster_stars, cluster_mass, 
                        cluster_radius, gas_particles, gas_mass, gas_radius, 
                        impact_parameter, v_infinity, timestep, filename)
            else:
                print filename + " already exists, skipping simulation"



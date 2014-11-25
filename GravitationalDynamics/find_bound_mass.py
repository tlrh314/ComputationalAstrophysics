from amuse.lab import *
from matplotlib import pyplot

converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.parsec)
G_si = converter.to_si(nbody_system.G)

#TODO: function with tidal radii?

#uses "HOP: A New Group Finding Algorithm for N-body Simulations"
def find_bound_bodies(bodies):
    return bodies.bound_subset(unit_converter=converter, G=G_si)

#assume that bodies that have higher kinetic than potential energy are
#not gravitationally bound, and returns bodies that are bound
def find_bound_bodies_using_energies(bodies):
    bound_particles = Particles()
    for body in bodies:
        body_energy = body.specific_kinetic_energy() + body.potential()
        if body_energy < 0 | body_energy.unit:
            bound_particles.add_particle(body)
    return bound_particles

if __name__ in '__main__':
    filename = "nbody.hdf5"
    bodies = read_set_from_file(filename, format="hdf5", copy_history=True)
    timestamps = []
    bound_masses = []
    for bodies_at_timestep in bodies.history:
        bound_bodies = find_bound_bodies(bodies_at_timestep)
        timestamp = bodies_at_timestep.get_timestamp()
        bound_mass = bound_bodies.total_mass()
        print "at timestamp", timestamp, "total mass of", bound_mass
        timestamps.append(timestamp)
        bound_masses.append(bound_mass)
    pyplot.scatter([i.number for i in timestamps], [i.number for i in bound_masses])
    pyplot.xlabel("time in MYear")
    pyplot.ylabel("bound mass in kg")
    pyplot.savefig("bound_mass_over_time")

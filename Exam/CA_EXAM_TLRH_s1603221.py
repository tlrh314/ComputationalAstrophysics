"""
Computational Astrophysics Exam A
Timo Halbesma, s1603221
January 14th, 2015. Version 2.0 (POST EXAM)

NB I still have to implement the rest of the plot function. In order to do so I
have to save the relative energy error as a function of time to file in the
hybrid_self_gravitating_cluster method. Then I can read the data in the
make_plots method.

Moreover, I have to implement the plot function for the
half-mass and core-radii as function of time.

Next, I have to rewrite the main method to run three runs and time those runs
in order to generate the table with the wall-clock times of the three runs.

Finally, I have to style my code a littlebit to make it more handsome.

Also, I still have to answer question 2.

Answers to questions.
Question 1) We prefer a split on mass because the Salpeter IMF is expected to
yield a mean mass of 0.333 MSun. This will result in assigning a higher number
of stars to the tree-code, and a lower number of stars to the direct N-Body
code. Since the direct N-Body code is O(N**2), a lower N is profitable. The
tree-code is O(n log n), thus, a higher N is managable in terms om computation
(wall clock) time. We do not choose to split on distance from the center of
mass because we want benefit from the higher precision of the direct N-Body
code in the entire Plummer sphere, not just inside/outside the center. We do
not split as a function of local density for the same reason.

Question 2)
"""


from time import time

from numpy import isnan

from amuse.units import units, nbody_system
from amuse.units.optparse import OptionParser

from amuse.support.console import set_printing_strategy
from amuse.io import write_set_to_file, read_set_from_file
from amuse.datamodel.particles import Channels

from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution

from amuse.couple.bridge import Bridge
from amuse.community.bhtree.interface import BHTree
from amuse.community.hermite0.interface import Hermite


def create_cluster(Ncl, Mstar, Rcl, Mcut):
    """
    Create a cluster with Salpeter IMF from 0.1 - Mstar | units.MSun,
    distributed in a virialzed Plummer sphere. Split stars into
    two particle sets: below_cut and above_cut for m < Mcut; m >= Mcut.
    """

    print "Deriving a set of", str(Ncl), "random masses",\
        "following a Salpeter IMF between 0.1 and",\
        str(Mstar), "MSun (alpha = -2.35)."
    salpeter_masses = new_salpeter_mass_distribution(Ncl, mass_max=Mstar)
    total_mass = salpeter_masses.sum()

    print "Distributing particles in a virialized Plummer sphere."
    nbody_converter = nbody_system.nbody_to_si(total_mass, Rcl)
    particles = new_plummer_model(Ncl, nbody_converter)
    particles.move_to_center()

    print "Setting masses of the stars in Plummer sphere."
    particles.mass = salpeter_masses

    print "Cluster virial radius", particles.virial_radius()

    # NB ParticlesSubset does not have the savepoint attribute?!
    # Without ,copy(), we have amuse.datamodel.particles.ParticlesSubset
    # With .copy(), we have amuse.datamodel.particles.Particles; save works (:
    stars_below_cut = particles.select(lambda m: m < Mcut, ["mass"])
    stars_above_cut = particles.select(lambda m: m >= Mcut, ["mass"])

    print "Number of stars below", str(Mcut), "=", len(stars_below_cut)
    print "Number of stars above", str(Mcut), "=", len(stars_above_cut)

    try:
        import Gnuplot
    except ImportError:
        print "Unable to produce plot: couldn't find Gnuplot."
    else:
        print "Plotting distribution of stars below and above cut-off."
        plotter = Gnuplot.Gnuplot()
        plotter('set terminal post enhance color')
        plotter('set output "CA_Exam_TLRH_PlummerSphere_SalpeterIMF_splitup.ps"')
        plotter.splot(stars_below_cut.position.value_in(units.parsec),
                      stars_above_cut.position.value_in(units.parsec),)

    # I return the particles superset because the particles superset is saved
    # when the subset is changed. This superset can be written to file. (:
    return particles, stars_below_cut, stars_above_cut, nbody_converter


def setup_codes(stars_below_cut, stars_above_cut, nbody_converter):
    if len(stars_below_cut) is 0:  # direct
        print "There are zero stars below the cut-off, using direct solver!"
        gravity = Hermite(nbody_converter)
        gravity.particles.add_particles(stars_above_cut)
        channels = Channels()
        channels.add_channel(gravity.particles.new_channel_to(stars_above_cut))
        return gravity, channels  # NB no bridge needed

    if len(stars_above_cut) is 0:  # tree
        print "There are zero stars above the cut-off, using tree solver!"
        gravity = BHTree(nbody_converter)
        gravity.particles.add_particles(stars_below_cut)
        channels = Channels()
        channels.add_channel(gravity.particles.new_channel_to(stars_below_cut))
        return gravity, channels  # NB no bridge needed

    gravity_low_mass = BHTree(nbody_converter)
    gravity_low_mass.particles.add_particles(stars_below_cut)

    gravity_high_mass = Hermite(nbody_converter, number_of_workers=4)
    gravity_high_mass.particles.add_particles(stars_above_cut)

    bridge = Bridge(timestep=0.01|units.Myr, use_threading=False)
    bridge.add_system(gravity_low_mass, (gravity_high_mass,))
    bridge.add_system(gravity_high_mass, (gravity_low_mass,))

    channels = Channels()
    channels.add_channel(gravity_low_mass.particles.new_channel_to(stars_below_cut))
    channels.add_channel(gravity_high_mass.particles.new_channel_to(stars_above_cut))

    return bridge, channels


def hybrid_self_gravitating_cluster(Ncl, Mstar, Rcl, t_end, dt,
                                    Mcut, filename):
    """ Create and evolve hybrid cluster; store data to file """

    stars_superset, stars_below_cut, stars_above_cut, nbody_converter =\
        create_cluster(Ncl, Mstar, Rcl, Mcut)

    bridge, channels = setup_codes(stars_below_cut, stars_above_cut,
                                   nbody_converter)

    time = 0. | units.Myr
    write_set_to_file(stars_superset.savepoint(time), filename, 'amuse',
            append_to_file=False)

    print "Saving data to", filename
    while time < (t_end - dt/2.):
        time += dt

        print "Evolving to", time
        # NB this is the solver, not the bridge when only one solver is used!
        # evolve_model, however, can be called both on bridge and solver (:
        bridge.evolve_model(time)
        channels.copy()

        write_set_to_file(stars_superset.savepoint(time),
            filename, 'amuse')


def update_progress(progress):
    """
    https://stackoverflow.com/questions/3160699/python-progress-bar

    update_progress() : Displays or updates a console progress bar

    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%
    """

    import sys

    # Modify length to change the length of the progress bar
    length = 42
    status = ""

    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt ):\r\n"
    if progress >= 1 - 0.0001:
        progress = 1
        status = "Done (:\r\n"

    block = int(round(length * progress))
    text = "\rprogress: [{0}] {1:.2f}% {2}".\
        format("#" * block + "-" * (length - block), progress * 100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def calculate_energy_error(particles, Etot_init):
    """ Given amuse.datamodel.particles.Particles & Etot_init, calculate dE """

    Ekin = particles.kinetic_energy()
    Epot = particles.potential_energy()
    Etot = Ekin + Epot

    dE = (Etot_init - Etot) / Etot  # loses unit, breaks plot

    if isnan(dE):
        dE = 0

    return dE


def make_plots(Ncl, Rcl, t_end):
    """ Not finished """
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot
    except ImportError:
        print "Unable to produce plots: couldn't find matplotlib"
    else:
        fig, ax = pyplot.subplots()
        dE, dt = [], []
        for data in ["hybrid", "tree", "direct"]:
            print "Generating plot data of {0} run".format(data)
            data_file = "CA_Exam_TLRH_{0}.amuse".format(data)
            stars_superset = read_set_from_file(data_file, 'amuse')

            # lim = abs(next(stars_superset.history).position).max()

            # Calculate initial energy error.
            Ekin_init = next(stars_superset.history).kinetic_energy()
            Epot_init = next(stars_superset.history).potential_energy()
            Etot_init = Ekin_init + Epot_init

            dE[:], dt[:] = [], []
            for snapshot in stars_superset.history:
                time = snapshot.get_timestamp().value_in(units.Myr)
                dt.append(time)
                dE.append(abs(calculate_energy_error(snapshot, Etot_init)))
                update_progress(float(time) / t_end.value_in(units.Myr))

            if data == "hybrid":
                pyplot.plot(dt, dE, 'b--', label='Hybrid')
            if data == "direct":
                pyplot.plot(dt, dE, 'ro', label='Direct')
            if data == "tree":
                pyplot.plot(dt, dE, 'g-', label='Tree')

        ax.set_xlabel("Time [Myr]")
        ax.set_ylabel("Relative Energy Error")
        ax.set_yscale('log', basex=2)
        ax.legend(loc=2)  # upper left
        ax.set_title(r'Cluster with $N=${0}, $r=${1}'.format(Ncl, Rcl) +
                     r', evolved until $t_{\rm end}$'+r'={0}'.format(t_end))
        pyplot.legend()
        pyplot.savefig("out.png")
        pyplot.close()


def parse_options():
    parser = OptionParser()
    parser.add_option("-N", dest="Ncl", type="int",
        default = 1000, help="number of stars [%default]")
    parser.add_option("-M", dest="Mstar", unit=units.MSun,
        type="float", default = 100|units.MSun,
        help="Maximum stellar mass [%default]")
    parser.add_option("-m", dest="Mcut", unit=units.MSun,
        type="float", default = 0.3|units.MSun,
        help="Cut-off mass [%default]")
    parser.add_option("-r", dest="Rcl", unit=units.parsec,
        type="float", default = 3|units.parsec,
        help="cluster virial radius [%default]")
    parser.add_option("-t", dest="dt", unit=units.Myr,
        type="float", default = 0.1|units.Myr,
        help="Time resolution of the simulation [%default]")
    parser.add_option("-T", dest="t_end", unit=units.Myr,
        type="float", default = 10|units.Myr,
        help="End time of the simulation [%default]")
    parser.add_option("-f", dest="filename",
        type="string", default="CA_Exam_TLRH_hybrid.amuse",
        help="The file in which to save the data [%default]")

    options, arguments = parser.parse_args()
    return options.__dict__


if __name__ == "__main__":
    set_printing_strategy("custom",
        preferred_units=[units.MSun, units.parsec, units.Myr, units.kms])

    options = parse_options()

    generate_data = False

    if generate_data:
        # Hybrid
        t_start = time()
        hybrid_self_gravitating_cluster(**options)
        print "Runtime of hybrid code =", time() - t_start

        # All stars above_cut, thus, direct
        t_start = time()
        options["Mcut"] = 0. | units.MSun
        options["filename"] = "CA_Exam_TLRH_direct.amuse"
        hybrid_self_gravitating_cluster(**options)
        print "Runtime of direct code =", time() - t_start

        # All stars below_cut, thus, tree
        t_start = time()
        options["Mcut"] = options["Mstar"]
        options["filename"] = "CA_Exam_TLRH_tree.amuse"
        hybrid_self_gravitating_cluster(**options) # Hybrid
        print "Runtime of tree code =", time() - t_start

    make_plots(options["Ncl"], options["Rcl"], options["t_end"])

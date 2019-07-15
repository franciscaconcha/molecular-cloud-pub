import numpy
from matplotlib import pyplot

from amuse.lab import *


def kroupa2001(masses):
    fit = []

    for mass in masses:
        if 0.01 <= mass < 0.08:
            fit.append(12.5 * numpy.power(mass, 0.7))
        elif 0.08 <= mass < 0.50:
            fit.append(numpy.power(mass, -0.3))
        elif 0.50 <= mass < 1.0:
            fit.append(0.5 * numpy.power(mass, -1.3))
        elif mass >= 1.0:
            fit.append(0.5 * numpy.power(mass, -1.3))

    print len(fit)

    return fit


def fill_mass_function_with_sink_mass(total_mass):
    print "Make mass function for M=", total_mass.in_(units.MSun)
    masses = [] | units.MSun

    while total_mass > 0 | units.MSun:
        #print "Available mass: ", total_mass.in_(units.MSun)
        mass = new_kroupa_mass_distribution(1, total_mass)[0]
        #print "Created 1 star of mass: ", mass.in_(units.MSun)
        if mass > total_mass:
            #print "Mass > available mass!"
            #print "Create a new star with the remaining available mass"
            mass = total_mass
        total_mass -= mass
        masses.append(mass)
        #print "New stars: ", len(masses)
        #print "*"

    print "N new stars=", len(masses), "total mass=", masses.sum().in_(units.MSun)
    return masses


def main():
    pyplot.style.use('paper')

    masses = numpy.arange(0.01, 100, 0.01)
    print len(masses)
    #fit = kroupa2001(masses)
    #pyplot.loglog(masses, fit, lw=3, c='red', label='Kroupa 2001')

    sink_masses = [1, 10, 100, 1000, 10000] | units.MSun
    colors = ['blue', 'lavender', 'pink', 'yellow', 'green']

    i = 0

    for s in sink_masses:
        for j in range(50):
            print "j = ", j
            sink_stars_masses = fill_mass_function_with_sink_mass(s).value_in(units.MSun)
            N = len(sink_stars_masses)
            fit = kroupa2001(numpy.sort(sink_stars_masses))
            kroupa_masses = new_kroupa_mass_distribution(N, mass_max=50 | units.MSun).value_in(units.MSun)  # Total mass is the mass of the sink
            print "total sink mass = {0} MSun, total kroupa = {1} MSun".format(sink_stars_masses.sum(),
                                                                               kroupa_masses.sum())

            #y, binEdges = numpy.histogram(kroupa_masses, bins=10)
            #bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
            pyplot.hist(sink_stars_masses, bins=100, color=colors[i], alpha=0.5)
            pyplot.loglog(numpy.sort(sink_stars_masses), fit, '-', c='red', lw=3, alpha=0.5)#, label="Kroupa 2001")

            #y, binEdges = numpy.histogram(sink_stars_masses, bins=10)
            #bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
            #pyplot.loglog(bincenters, y, '-', c=colors[i])#, lw=3, label="Kroupa 2001")

        i += 1
        pyplot.legend()

        pyplot.show()


if __name__ in ("__main__", "__plot__"):
    main()

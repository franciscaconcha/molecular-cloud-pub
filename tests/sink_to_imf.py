import numpy
from matplotlib import pyplot
from decimal import Decimal

from amuse.lab import *


def kroupa2001(masses):
    imf = []

    for mass in masses:
        if 0.01 <= mass < 0.08:
            imf.append(12.5 * numpy.power(mass, 0.7))
        elif 0.08 <= mass < 0.50:
            imf.append(numpy.power(mass, -0.3))
        elif 0.50 <= mass < 1.0:
            imf.append(0.5 * numpy.power(mass, -1.3))
        elif mass >= 1.0:
            imf.append(0.5 * numpy.power(mass, -1.3))

    return imf


def salpeter1955(masses):
    return numpy.power(masses, -2.35)


def fill_mass_function_with_sink_mass(total_mass, imf):
    #print "Make mass function for M=", total_mass.in_(units.MSun)
    masses = [] | units.MSun

    if imf == 'Kroupa':
        while total_mass > 0 | units.MSun:
            mass = new_kroupa_mass_distribution(1, mass_max=total_mass)[0]
            if mass > total_mass:
                    mass = total_mass
            total_mass -= mass
            masses.append(mass)

    elif imf == 'Salpeter':
        while total_mass > 0 | units.MSun:
            mass = new_salpeter_mass_distribution(1, mass_max=total_mass)[0]
            if mass > total_mass:
                mass = total_mass
            total_mass -= mass
            masses.append(mass)

    else:
        print "Error in IMF name"
        return 0

    #mean_mass = masses.sum() / len(masses)
    #Nstars = int(total_mass/mean_mass)
    #masses = new_kroupa_mass_distribution(Nstars, 100|units.MSun)
    #print "Generate N=", Nstars, "stars of mean mass m=", mean_mass.in_(units.MSun)
    #print "total mass generated M=", masses.sum().in_(units.MSun)

    #print "N new stars=", len(masses), "total mass=", masses.sum().in_(units.MSun)
    return masses


def main(imf):
    pyplot.style.use('paper')
    fig = pyplot.figure()

    sink_masses = [1, 10, 100, 1000, 10000] | units.MSun

    for s in sink_masses:
        pyplot.clf()
        Ns = []
        for j in range(1000):
            sink_stars_masses = fill_mass_function_with_sink_mass(s, imf).value_in(units.MSun)
            sink_stars_masses = sink_stars_masses[sink_stars_masses >= 0.01]  # Lower limit for Kroupa IMF
            Ns.append(len(sink_stars_masses))

            if imf == 'Kroupa':
                fit = kroupa2001(numpy.sort(sink_stars_masses))
            elif imf == 'Salpeter':
                fit = salpeter1955(numpy.sort(sink_stars_masses))
            else:
                print "Error in IMF name"
                return 0

            pyplot.loglog(numpy.sort(sink_stars_masses), fit, '-', c='black', lw=3, alpha=0.2)

        if imf == 'Kroupa':
            kroupa_masses = new_kroupa_mass_distribution(int(numpy.mean(Ns)), mass_max=s).value_in(units.MSun)
            kroupa_imf = kroupa2001(numpy.sort(kroupa_masses))
            pyplot.loglog(numpy.sort(kroupa_masses), kroupa_imf, c='red', lw=3)
        elif imf == 'Salpeter':
            salpeter_masses = new_salpeter_mass_distribution(int(numpy.mean(Ns)), mass_max=s).value_in(units.MSun)
            salpeter_imf = salpeter1955(numpy.sort(salpeter_masses))
            pyplot.loglog(numpy.sort(salpeter_masses), salpeter_imf, c='red', lw=3)
        else:
            print "Error in IMF name"
            return 0

        pyplot.xlabel(r'$m$ [$\mathrm{M}_{\odot}$]')
        pyplot.ylabel(r'$\xi(m)\ dm$')
        plot_title = '{0} IMF, '.format(imf) + \
                     r'sink mass {0}'.format(Decimal(s.value_in(units.MSun))) + \
                     r'$\mathrm{M}_{\odot}$, $<N_*>$ = ' + '{0}'.format(int(numpy.mean(Ns)))
        pyplot.title(plot_title)
        pyplot.legend()
        #pyplot.show()
        pyplot.savefig('tests/results/IMF/{0}_sink_{1}MSun.png'.format(imf, s.value_in(units.MSun)))




def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-i", dest="imf",
                      default="Kroupa",
                      help="IMF [Kroupa, Salpeter]")

    return result

if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)

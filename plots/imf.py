import numpy
from matplotlib import pyplot
from matplotlib.lines import Line2D
import os

from amuse.lab import *

from mycolors import *

def imf(path, save_path, N, nruns=1):
    fig = pyplot.figure()
    ax = pyplot.gca()

    # Create a Kroupa 2001 distribution with the same number of stars, to compare
    IMF_masses = new_kroupa_mass_distribution(N,
                                              mass_max=150.0 | units.MSun).value_in(units.MSun)

    #kroupa_limits = [0.01, 0.08, 0.5, 1.0, 1.9, 150.0]  # Added the 1.9 MSun 'bin' for photoevap limit
    kroupa_limits = [0.08, 0.5, 1.0, 1.9, 150.0]  # Added the 1.9 MSun 'bin' for photoevap limit

    for r in range(nruns):
        filepath = '{0}/{1}/gravity_stars.hdf5'.format(path,
                                                       r)

        stars = read_set_from_file(filepath, "hdf5", close_file=True)
        star_masses = stars.stellar_mass.value_in(units.MSun)

        ax.hist(star_masses, bins=kroupa_limits, histtype=u'step', edgecolor='k', lw=3, alpha=0.5)

    ax.hist(IMF_masses, bins=kroupa_limits, histtype=u'step', edgecolor='r', lw=3)

    lines = [Line2D([0], [0], color='red', linewidth=3),
             Line2D([0], [0], color='k', linewidth=3, alpha=0.5)]
    labels = ['Kroupa IMF', 'Simulations']
    pyplot.legend(lines, labels)

    # 1.9 MSun limit, for photoevaporation
    ax.axvline(1.9,
               lw=3,
               ls='--',
               c='blue')
    ax.text(1.5,
            3000,
            r'$M_* = 1.9 M_{{\odot}}$',
            color='blue',
            #transform=ax.transAxes,
            rotation=90)

    ax.set_xlabel(r'$M_*$ [$\mathrm{M}_{\odot}$]')
    ax.set_ylabel(r'$N_*$')
    #ax.set_title(r'Total $N_*$ = {0}'.format(N))

    #ax.set_xlim([0.001, 150])

    xtick_labels = [item.get_text() for item in pyplot.gca().get_xticklabels()]
    xtick_labels[0] = '0.08'

    pyplot.gca().set_xticklabels(xtick_labels)

    ax.set_xscale('log')

    pyplot.show()
    #pyplot.savefig('{0}/IMF_vs_simulation.png'.format(save_path))


def colortest():
    x = numpy.linspace(0, 2 * numpy.pi, 64)
    y = numpy.cos(x)
    y = numpy.cos(x)
    lab = ['0', '1', '2', '3', '4']
    fig = pyplot.figure()
    for i in range(len(lab)):
        pyplot.plot(x, i*y, c=colors[lab[i]], lw=5, label=lab[i])
    pyplot.legend()
    pyplot.show()


def main(path, save_path, tend, dt_diag, Ncloud, Mcloud, Rcloud):
    # My own style sheet, comment out if not needed
    pyplot.style.use('paper')

    #colortest()

    imf(path, save_path, 10000, nruns=12)



def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-p", dest="path",
                      default=".",
                      help="input filename")
    result.add_option("-s", dest="save_path",
                      default="./results/",
                      help="save path for results")
    result.add_option("--tend", dest="tend",
                      unit=units.Myr,
                      type="float",
                      default=5.0 | units.Myr,
                      help="end time")
    result.add_option("--dt_diag", dest="dt_diag",
                      unit=units.Myr,
                      type="float",
                      default=0.1 | units.Myr,
                      help="diagnosticstime step")
    result.add_option("--Ncloud", dest="Ncloud",
                      default=4000,
                      type="float",
                      help="number of gas particles.")
    result.add_option("--Mcloud", dest="Mcloud",
                      unit=units.MSun,
                      type="float",
                      default=10000 | units.MSun,
                      help="cloud mass")
    result.add_option("--Rcloud", dest="Rcloud",
                      unit=units.parsec,
                      type="float",
                      default=2 | units.parsec,
                      help="cloud size")

    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    numpy.random.seed(3141)
    main(**o.__dict__)

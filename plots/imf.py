import numpy
from matplotlib import pyplot
from matplotlib.lines import Line2D
import os

from amuse.lab import *

from mycolors import *


def imf(path, save_path, N, nruns=1, save=False):
    fig = pyplot.figure()
    ax = pyplot.gca()

    # Create a Kroupa 2001 distribution with the same number of stars, to compare
    IMF_masses = new_kroupa_mass_distribution(N,
                                              mass_max=150.0 | units.MSun).value_in(units.MSun)

    kroupa_limits = [0.08, 0.5, 1.0, 1.9, 150.0]  # Added the 1.9 MSun 'bin' for photoevap limit

    all_mean_masses = []

    for r in range(nruns):
        filepath = '{0}/{1}/gravity_stars.hdf5'.format(path,
                                                       r)

        stars = read_set_from_file(filepath, "hdf5", close_file=True)
        star_masses = stars.stellar_mass.value_in(units.MSun)
        all_mean_masses.append(numpy.mean(star_masses))

        ax.hist(star_masses, bins=kroupa_limits, histtype=u'step', edgecolor='gray', lw=3, alpha=0.5)

    ax.hist(IMF_masses, bins=kroupa_limits, histtype=u'step', edgecolor=colors['red'], lw=4)

    print numpy.mean(all_mean_masses), numpy.std(all_mean_masses)

    lines = [Line2D([0], [0], color=colors['red'], linewidth=3),
             Line2D([0], [0], color='gray', linewidth=3, alpha=0.5)]
    labels = ['Kroupa IMF', 'Simulations']
    pyplot.legend(lines, labels)

    # 1.9 MSun limit, for photoevaporation
    ax.axvline(1.9,
               lw=3,
               ls='--',
               c=colors['blue'])
    ax.text(1.5,
            3000,
            r'$M_* = 1.9 M_{{\odot}}$',
            color=colors['blue'],
            #transform=ax.transAxes,
            rotation=90)

    ax.set_xlabel(r'$\mathrm{M}_*$ [$\mathrm{M}_{\odot}$]')
    ax.set_ylabel(r'$\mathrm{N}_*$')
    #ax.set_title(r'Total $N_*$ = {0}'.format(N))

    #ax.set_xlim([0.001, 150])

    xtick_labels = [item.get_text() for item in pyplot.gca().get_xticklabels()]
    xtick_labels[0] = '0.08'

    pyplot.gca().set_xticklabels(xtick_labels)

    ax.set_xscale('log')

    if save:
        pyplot.savefig('{0}/IMF_vs_simulation.png'.format(save_path))
    else:
        pyplot.show()


def colortest():
    x = numpy.linspace(0, 2 * numpy.pi, 64)
    y = numpy.cos(x)
    lab = [str(a) for a in range(10)]
    fig = pyplot.figure()
    for i in range(len(lab)):
        pyplot.plot(x, i * y, c=runcolors[int(lab[i])], lw=5, label=lab[i])
    pyplot.legend()
    pyplot.show()


def main(path, save_path, nruns, save):
    # My own style sheet, comment out if not needed
    pyplot.style.use('paper')

    #colortest()

    imf(path, save_path, 10000, nruns, save)


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-p", dest="path",
                      default=".",
                      help="input filename")
    result.add_option("-s", dest="save_path",
                      default="./results/",
                      help="save path for results")
    result.add_option("-n", dest="nruns",
                      type="int",
                      default=10,
                      help="number of runs to plot")
    result.add_option("-S", dest="save",
                      type="int",
                      default=0,
                      help="if 1, save figure")

    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    numpy.random.seed(3141)
    main(**o.__dict__)

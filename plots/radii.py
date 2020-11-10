import numpy
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
import os

from amuse.lab import *
from amuse import io

from mycolors import *
from legends import *


def Rvir(open_path, save_path, nruns, save):
    """ Virial radius in time.

    :param open_path: path to results file
    :param N: number of stars in results
    :param save_path: path to save figure
    :param t_end: final time to use when plotting
    :param save: if True, figure will be saved
    :param nrun: run number to use for the plot
    """
    fig = pyplot.figure()

    times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    Rvir = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

    for n in range(nruns):
        path = '{0}/{1}/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if 'hydro_stars' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        for f in files[1:-1]:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            #t = float(f.split('t')[1].split('.hdf5')[0])
            t = stars.get_timestamp().value_in(units.Myr)
            times[n].append(t)

            vr = stars.virial_radius().value_in(units.parsec)
            Rvir[n].append(vr)

    for n in range(nruns):
        pyplot.plot(times[n], Rvir[n], c=runcolors[n], lw=3, label='Run {0}'.format(n))

    pyplot.legend(loc='best', ncol=2)
    pyplot.xlabel('Time [Myr]')
    pyplot.ylabel(r'$\mathrm{R}_\mathrm{vir}$ [pc]')

    if save:
        pyplot.savefig('{0}/Rvir.png'.format(save_path))
    else:
        pyplot.show()


def Rhm(open_path, save_path, nruns, save):
    fig = pyplot.figure()

    times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    Rhm = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

    for n in range(nruns):
        path = '{0}/{1}/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if 'hydro_stars' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        for f in files[1:-1]:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            # t = float(f.split('t')[1].split('.hdf5')[0])
            t = stars.get_timestamp().value_in(units.Myr)

            converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), 3.0 | units.parsec)

            r_hm, mf = stars.LagrangianRadii(mf=[0.5], unit_converter=converter)
            hmr = r_hm.value_in(units.parsec)[0]

            times[n].append(t)
            Rhm[n].append(hmr)

    for n in range(nruns):
        pyplot.plot(times[n], Rhm[n], c=runcolors[n], lw=3, label='Run {0}'.format(n))

    pyplot.legend(loc='best', ncol=2)

    pyplot.xlabel('Time [Myr]')

    pyplot.ylabel(r'$\mathrm{R}_\mathrm{hm}$ [pc]')

    if save:

        pyplot.savefig('{0}/Rhm.png'.format(save_path))

    else:

        pyplot.show()


def main(open_path, N, save_path, t_end, save, distance, nruns, time):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    Rvir(open_path, save_path, nruns, save)
    Rhm(open_path, save_path, nruns, save)


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()

    result.add_option("-p", dest="open_path", type="string", default='/media/fran/data1/photoevap/results',
                      help="path to results to plot [%default]")
    result.add_option("-N", dest="N", type="int", default=1000,
                      help="number of stars [%default]")
    result.add_option("-S", dest="save", type="int", default=0,
                      help="save plot? [%default]")
    result.add_option("-s", dest="save_path", type="string", default='/media/fran/data1/photoevap-paper/figures',
                      help="path to save the results [%default]")
    result.add_option("-e", dest="t_end", type="float", default='5.0',
                      help="end time to use for plots [%default]")
    result.add_option("-d", dest="distance", type="float", default=0.0,
                      help="When using galactic potential, ('projected') distance to galactic center [%default]")
    result.add_option("-n", dest="nruns", type="int", default=1,
                      help="number of runs to plot for averages [%default]")
    result.add_option("-t", dest="time", type="int", default=0,
                      help="make movie? [%default]")
    return result


if __name__ == '__main__':
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)

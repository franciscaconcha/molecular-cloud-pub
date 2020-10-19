import numpy
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
import os

from amuse.lab import *
from amuse import io

from mycolors import *
from legends import *
from readradii import *


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

    min_times, max_times = [], []

    for n in range(nruns):
        path = '{0}/{1}/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))
        min_times.append(float(files[1].split('t')[1].split('.hdf5')[0]))  # [0] is gravity_stars.hdf5
        max_times.append(float(files[-1].split('t')[1].split('.hdf5')[0]))

    dt = 0.005
    times = numpy.arange(1.0,
                         max(max_times),
                         dt)

    print min(min_times), max(max_times)

    for n in range(nruns):
        N = Ns[n]
        virial_radii = readRvir[n]
        """virial_radii = []

        for t in times:
            f = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, n, N, t)
            try:
                stars = io.read_set_from_file(f, 'hdf5', close_file=True)

                born_stars = stars[stars.born]
                converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), 3.0 | units.parsec)

                vr = born_stars.virial_radius().value_in(units.parsec)
                virial_radii.append(vr)
                #try:
                #    r_hm, mf = born_stars.LagrangianRadii(mf=[0.5], unit_converter=converter)
                #    hmr = r_hm.value_in(units.parsec)[0]
                #except:
                #    hmr = 0.0
            except io.base.IoException:
                virial_radii.append(numpy.nan)

        print n
        print virial_radii"""
        pyplot.plot(times, virial_radii, c=runcolors[n], lw=3, label='Run {0}'.format(n))

    pyplot.legend(loc='best', ncol=2)
    pyplot.xlabel('Time [Myr]')
    pyplot.ylabel(r'$\mathrm{R}_\mathrm{vir}$ [pc]')

    if save:
        pyplot.savefig('{0}/Rvir.png'.format(save_path))
    else:
        pyplot.show()


def Rhm(open_path, save_path, nruns, save):
    fig = pyplot.figure()

    min_times, max_times = [], []

    for n in range(nruns):
        path = '{0}/{1}/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))
        min_times.append(float(files[1].split('t')[1].split('.hdf5')[0]))  # [0] is gravity_stars.hdf5
        max_times.append(float(files[-1].split('t')[1].split('.hdf5')[0]))

    dt = 0.005
    times = numpy.arange(1.0,
                         max(max_times),
                         dt)

    print min(min_times), max(max_times)
    for n in range(nruns):
        N = Ns[n]
        hm_radii = readRhm[n]
        """hm_radii = []

        for t in times:
            f = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, n, N, t)
            try:
                stars = io.read_set_from_file(f, 'hdf5', close_file=True)

                born_stars = stars[stars.born]
                converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), 3.0 | units.parsec)

                r_hm, mf = born_stars.LagrangianRadii(mf=[0.5], unit_converter=converter)
                hmr = r_hm.value_in(units.parsec)[0]

                hm_radii.append(hmr)

            except io.base.IoException:
                hm_radii.append(numpy.nan)

        print n
        print hm_radii"""

        pyplot.plot(times, hm_radii, c=runcolors[n], lw=3, label='Run {0}'.format(n))

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

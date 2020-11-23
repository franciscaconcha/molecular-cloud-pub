import numpy
from matplotlib import pyplot
from scipy import stats
from sklearn.neighbors import KDTree
import os

from amuse.lab import *
from amuse import io

#movie command
#"ffmpeg -framerate 5 -i {0}/%01d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {0}/movie.mp4

from legends import *  # My own custom legend definitions
from mycolors import *


def get_mean(data):
    """Returns sorted mean of the data"""
    try:
        mean_data = numpy.mean(data, axis=0)
        #mean_cumulative = numpy.mean(cumulative, axis=0)
        devs = numpy.std(data, axis=0)
    except:
        max_len = 0
        for a in data:
            if len(a) > max_len:
                max_len = len(a)

        new_sorted = []
        for a in data:
            b = numpy.pad(a, (max_len - len(a), 0), 'constant',#)
                constant_values=(min([min(r) for r in data])))
            new_sorted.append(b)

        mean_data = numpy.mean(new_sorted, axis=0)
        devs = numpy.std(new_sorted, axis=0)

        max_len = 0
        """for a in cumulative:
            if len(a) > max_len:
                max_len = len(a)

        new_sorted = []
        for a in cumulative:
            b = numpy.pad(a, (max_len - len(a), 0), 'constant')
            # constant_values=(min([min(r) for r in all_initial])))
            new_sorted.append(b)

        mean_cumulative = numpy.mean(new_sorted, axis=0)
        devs = numpy.std(new_sorted, axis=0)"""

    mean_data = numpy.sort(mean_data)

    return mean_data, devs


def mean_mass_in_time(open_path, save_path, nruns, save):

    fig = pyplot.figure()

    times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    stdev = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

    for n in range(nruns):
        path = '{0}/{1}/disks/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if 'N' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        for f in files:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            t = float(f.split('t')[1].split('.hdf5')[0])
            times[n].append(t)

            stars = stars[stars.disked]
            disk_masses = stars.disk_mass.value_in(units.MJupiter)
            masses[n].append(numpy.mean(disk_masses))
            stdev[n].append(numpy.std(disk_masses))

    # Indices for time in which star formation ends for each run
    sf_end_indices = [85, 86, 224, 127, 196, 171]

    for n in range(nruns):
        i = sf_end_indices[n]
        pyplot.plot(times[n][:i],
                    masses[n][:i],
                    c=runcolors[n],
                    lw=3,
                    label="Run {0}".format(n))
        pyplot.plot(times[n][i:],
                    masses[n][i:],
                    c=runcolors[n],
                    lw=3,
                    ls=":",
                    )
        if n == 0 or n == 2:
            pyplot.fill_between(times[n],
                                numpy.array(masses[n]) - numpy.array(stdev[n]),
                                numpy.array(masses[n]) + numpy.array(stdev[n]),
                                facecolor=runcolors[n],
                                alpha=0.2,
                                )

    pyplot.legend()
    pyplot.xlabel(r'Time [Myr]')
    pyplot.ylabel(r'Mean disc mass [$\mathrm{M}_{Jup}$]')

    if save:
        pyplot.savefig('{0}/mean_mass_vs_time.png'.format(save_path))


def delta_vs_time(open_path, save_path, nruns, save):
    """
    delta = Mgas/Mdust
    :param open_path:
    :param save_path:
    :param nruns:
    :param save:
    :return:
    """

    fig = pyplot.figure()

    times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    deltas = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    stdev = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

    for n in range(nruns):
        path = '{0}/{1}/disks/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if 'N' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        for f in files:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            t = float(f.split('t')[1].split('.hdf5')[0])
            times[n].append(t)

            stars = stars[stars.disked]
            dust_masses = stars.disk_dust_mass.value_in(units.MJupiter)
            gas_masses = stars.disk_gas_mass.value_in(units.MJupiter)
            deltas[n].append(numpy.mean(dust_masses)/numpy.mean(gas_masses))
            stdev[n].append(numpy.std(dust_masses)/numpy.std(gas_masses))

    ds, ts = [], []
    for n in range(nruns):
        ds.append(deltas[n])
        ts.append(times[n])

    meants, devts = get_mean(ts)
    meands, devds = get_mean(ds)

    pyplot.plot(meants,
                meands)
    pyplot.fill_between(meants,
                        meands - devds,
                        meands + devds,
                        alpha=0.5)

    # Indices for time in which star formation ends for each run
    """sf_end_indices = [85, 86, 224, 127, 196, 171]

    for n in range(nruns):
        i = sf_end_indices[n]
        pyplot.plot(times[n][:i],
                    deltas[n][:i],
                    c=runcolors[n],
                    lw=3,
                    label="Run {0}".format(n))
        pyplot.plot(times[n][i:],
                    deltas[n][i:],
                    c=runcolors[n],
                    lw=3,
                    ls=":",
                    )"""
    """if n == 0 or n == 2:
            pyplot.fill_between(times[n],
                                numpy.array(deltas[n]) - numpy.array(stdev[n]),
                                numpy.array(deltas[n]) + numpy.array(stdev[n]),
                                facecolor=runcolors[n],
                                alpha=0.2,
                                )"""

    #pyplot.legend()

    #pyplot.yscale('log')

    pyplot.xlabel(r'Time [Myr]')
    pyplot.ylabel(r'$M_{dust} / M_{gas}$')

    if save:
        pyplot.savefig('{0}/delta_vs_time.png'.format(save_path))



def main(open_path, N, save_path, t_end, save, nruns):

    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    #mean_mass_in_time(open_path, save_path, nruns, save)
    delta_vs_time(open_path, save_path, nruns, save)

    if not save:
        pyplot.show()


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
    result.add_option("-n", dest="nruns", type="int", default=1,
                      help="number of runs to plot for averages [%default]")
    return result


if __name__ == '__main__':
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)

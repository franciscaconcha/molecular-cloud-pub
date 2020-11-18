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


def mean_mass_in_time(open_path, save_path, nruns):

    fig = pyplot.figure()

    times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

    for n in range(nruns):
        path = '{0}/{1}/disks/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        for f in files:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            t = float(f.split('t')[1].split('.hdf5')[0])
            times[n].append(t)

            # want to consider all stars to ever have a disk,
            # even if it has been dispersed (disked == False)
            always_disked = stars[stars.stellar_mass <= 1.9 | units.MSun]
            disk_masses = always_disked.disk_mass.value_in(units.MJupiter)
            masses[n].append(numpy.mean(disk_masses))

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

    pyplot.legend()
    pyplot.xlabel(r'Time [Myr]')
    pyplot.ylabel(r'Mean disc mass [$\mathrm{M}_{Jup}$]')


def main(open_path, N, save_path, t_end, save, nruns):

    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    mean_mass_in_time(open_path, save_path, nruns)

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

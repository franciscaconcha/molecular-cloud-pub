import numpy
from matplotlib import pyplot
import os

from amuse.lab import *
from amuse import io

# movie command
# "ffmpeg -framerate 5 -i {0}/%01d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {0}/movie.mp4


def stars(open_path, N, save_path, t_end, save, nrun):
    """ Function to create Figure 1 on the paper.

    :param open_path: path to results file
    :param N: number of stars in results
    :param save_path: path to save figure
    :param t_end: final time to use when plotting
    :param save: if True, figure will be saved
    :param nrun: run number to use for the plot
    """
    fig = pyplot.figure(figsize=(10, 10))
    ax = fig.gca()

    dt = 0.01
    times = numpy.arange(0.715, t_end + dt, dt)

    i = 0

    for t in times:
        f = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, nrun, N, t)
        stars = io.read_set_from_file(f, 'hdf5', close_file=True)
        print len(stars), len(stars[stars.dispersed])
        #stars = stars[stars.tborn.value_in(units.Myr) <= t]

        pyplot.clf()

        pyplot.scatter(stars.x.value_in(units.parsec),
                       stars.y.value_in(units.parsec),
                       marker='o',
                       s=stars.disk_radius.value_in(units.au),
                       c='gray',
                       alpha=0.5,
                       lw=1)

        pyplot.xlim([-2.5, 2.5])
        pyplot.ylim([-2.5, 2.5])

        pyplot.xlabel("x [pc]")
        pyplot.ylabel("y [pc]")

        pyplot.axes().set_aspect('equal')

        pyplot.savefig('{0}/{1}.png'.format(save_path, i))

        print "{0}, t = {1} Myr, N = {2}".format(i, t, len(stars))

        i += 1


def radius_vs_time(open_path, N, save_path, t_end, save, nrun):
    last_file = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, nrun, N, t_end)
    last_stars = io.read_set_from_file(last_file, 'hdf5', close_file=True)

    fig = pyplot.figure()
    ax = fig.gca()

    dt = 0.01

    times = numpy.arange(0.715, t_end, dt)

    files = os.listdir('{0}/{1}'.format(open_path, nrun))  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'

    sink_files = [x for x in files if '.hdf5' in x]

    sink_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    """times = []
    for sf in sink_files:
        #print sf.split('t')[1].split('.hdf5')
        t = sf.split('t')[1].split('.hdf5')[0]
        times.append(float(t))"""

    radius_dict, times_dict = {}, {}

    last_file = '{0}/{1}/{2}'.format(open_path, nrun, sink_files[-1])

    last_stars = io.read_set_from_file(last_file, 'hdf5', close_file=True)
    #last_stars = last_stars[last_stars.tborn.value_in(units.Myr) <= t_end]

    #print last_stars.tborn.in_(units.Myr)

    for k in last_stars.key:
        radius_dict[k] = []
        times_dict[k] = []

    for t in times[1:]:
        print t
        f = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, nrun, N, t)
        stars = io.read_set_from_file(f, 'hdf5', close_file=True)
        stars = stars[stars.tborn.value_in(units.Myr) <= t]
        stars = stars[stars.disked]
        #print len(stars), len(stars[stars.dispersed])
        #print stars.disk_radius.in_(units.au)
        for s in stars:
            #radius_dict[s.key].append(s.disk_radius.value_in(units.au))
            radius_dict[s.key].append(s.disk_mass.value_in(units.MJupiter))
            times_dict[s.key].append(t)

    for k in last_stars.key:
        pyplot.plot(times_dict[k], radius_dict[k])

    pyplot.show()


def main(open_path, N, save_path, t_end, save, Rvir, distance, nruns, movie):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    #stars(open_path, N, save_path, t_end, save, nruns)
    radius_vs_time(open_path, N, save_path, t_end, save, nruns)


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
    result.add_option("-R", dest="Rvir", type="float",
                      unit=units.parsec, default=0.5,
                      help="cluster virial radius [%default]")
    result.add_option("-d", dest="distance", type="float", default=0.0,
                      help="When using galactic potential, ('projected') distance to galactic center [%default]")
    result.add_option("-n", dest="nruns", type="int", default=1,
                      help="number of runs to plot for averages [%default]")
    result.add_option("-m", dest="movie", type="int", default=0,
                      help="make movie? [%default]")
    return result


if __name__ == '__main__':
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)

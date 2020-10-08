import numpy
from matplotlib import pyplot
import matplotlib.gridspec as gridspec

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
    fig = pyplot.figure()
    ax = pyplot.gca()

    f = '{0}/{1}/gravity_stars.hdf5'.format(open_path, nrun)
    stars = io.read_set_from_file(f, 'hdf5', close_file=True)
    stars.mass = stars.stellar_mass
    print len(stars)

    disked_stars = stars[stars.stellar_mass <= 1.9 | units.MSun]
    massive_stars = stars[stars.stellar_mass > 1.9 | units.MSun]

    ax.scatter(disked_stars.x.value_in(units.parsec),
               disked_stars.y.value_in(units.parsec),
               marker='o',
               #s=stars.disk_radius.value_in(units.au),
               c='gray',
               alpha=0.5,
               lw=1)

    ax.scatter(massive_stars.x.value_in(units.parsec),
               massive_stars.y.value_in(units.parsec),
               marker='*',
               #s=stars.disk_radius.value_in(units.au),
               c='red',
               alpha=0.5,
               lw=1)

    ax.set_xlabel("x [pc]")
    ax.set_ylabel("y [pc]")

    ax.set_aspect('equal')

    fig.suptitle("Run \#{0}, N={1}".format(nrun, len(stars)))
    pyplot.show()


def all_runs(open_path, N, save_path, t_end, save, nruns):
    fig, axs = pyplot.subplots(4, 3, subplot_kw=dict(aspect='equal', adjustable='box-forced'))  # rows, columns

    ax = {0: axs[0, 0],
          1: axs[0, 1],
          2: axs[0, 2],
          3: axs[1, 0],
          4: axs[1, 1],
          5: axs[1, 2],
          6: axs[2, 0],
          7: axs[2, 1],
          8: axs[2, 2],
          9: axs[3, 0],
          10: axs[3, 1],
          11: axs[3, 2]
    }

    for n in range(nruns):
        f = '{0}/{1}/gravity_stars.hdf5'.format(open_path, n)
        stars = io.read_set_from_file(f, 'hdf5', close_file=True)

        disked_stars = stars[stars.stellar_mass <= 1.9 | units.MSun]
        massive_stars = stars[stars.stellar_mass > 1.9 | units.MSun]

        print n + 1, len(stars), stars.Qparameter(), stars.box_counting_dimension()

        """ax[n].scatter(disked_stars.x.value_in(units.parsec),
                       disked_stars.y.value_in(units.parsec),
                       marker='o',
                       #s=stars.disk_radius.value_in(units.au),
                       c='gray',
                       alpha=0.5,
                       lw=1)

        ax[n].scatter(massive_stars.x.value_in(units.parsec),
                       massive_stars.y.value_in(units.parsec),
                       marker='*',
                       #s=stars.disk_radius.value_in(units.au),
                       c='red',
                       alpha=0.5,
                       lw=1)

        #ax0.set_xlabel("x [pc]")
        #ax0.set_ylabel("y [pc]")

        #ax[n].set_xlim([-3, 3])
        #ax[n].set_ylim([-3, 3])

    pyplot.show()"""


def time_ev(open_path, N, save_path, t_end, save, nruns):
    """ Function to create Figure 1 on the paper.

    :param open_path: path to results file
    :param N: number of stars in results
    :param save_path: path to save figure
    :param t_end: final time to use when plotting
    :param save: if True, figure will be saved
    :param nrun: run number to use for the plot
    """
    #ax0 = pyplot.subplot(311)
    #ax1 = pyplot.subplot(312)
    #ax2 = pyplot.subplot(313)

    for n in range(nruns):
        pyplot.clf()
        fig = pyplot.figure(figsize=(15, 10))
        grid = gridspec.GridSpec(ncols=4, nrows=2, figure=fig)
        ax0 = fig.add_subplot(grid[:2, :2])  # row, column
        ax1 = fig.add_subplot(grid[-2, 2:])
        ax2 = fig.add_subplot(grid[-1, 2:])

        f = '{0}/{1}/gravity_stars.hdf5'.format(open_path, n)
        stars = io.read_set_from_file(f, 'hdf5', close_file=True)
        stars.mass = stars.stellar_mass
        print len(stars)

        dt = 0.05

        tmin = 20 | units.Myr
        tmax = 0.05 | units.Myr
        for s in stars:
            if s.tborn < tmin:
                tmin = s.tborn
            elif s.tborn > tmax:
                tmax = s.tborn

        print "first stars at ", tmin.in_(units.Myr)
        print "last stars at ", tmax.in_(units.Myr)

        times = numpy.arange(tmin.value_in(units.Myr),
                             tmax.value_in(units.Myr),
                             dt)

        disked_stars = stars[stars.stellar_mass <= 1.9 | units.MSun]
        massive_stars = stars[stars.stellar_mass > 1.9 | units.MSun]

        ax0.scatter(disked_stars.x.value_in(units.parsec),
                       disked_stars.y.value_in(units.parsec),
                       marker='o',
                       #s=stars.disk_radius.value_in(units.au),
                       c='gray',
                       alpha=0.5,
                       lw=1)

        ax0.scatter(massive_stars.x.value_in(units.parsec),
                       massive_stars.y.value_in(units.parsec),
                       marker='*',
                       #s=stars.disk_radius.value_in(units.au),
                       c='red',
                       alpha=0.5,
                       lw=1)

        ax0.set_xlabel("x [pc]")
        ax0.set_ylabel("y [pc]")

        ax0.set_aspect('equal')

        virial_radii, hm_radii = [], []

        for t in times:
            born_stars = stars[stars.tborn.value_in(units.Myr) <= t]
            converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), 3.0 | units.parsec)

            vr = born_stars.virial_radius().value_in(units.parsec)

            try:
                r_hm, mf = born_stars.LagrangianRadii(mf=[0.5], unit_converter=converter)
                hmr = r_hm.value_in(units.parsec)[0]
            except:
                hmr = 0.0

            virial_radii.append(vr)
            hm_radii.append(hmr)

        #ax1.plot(times, virial_radii)
        ax2.plot(times, hm_radii)

        #pyplot.xlim([-2.5, 2.5])
        #pyplot.ylim([-2.5, 2.5])

        #pyplot.xlabel("x [pc]")
        #pyplot.ylabel("y [pc]")

        #pyplot.axes().set_aspect('equal')

        #pyplot.savefig('{0}/{1}.png'.format(save_path, i))

        #print "{0}, t = {1} Myr, N = {2}".format(i, t, len(stars))
        fig.suptitle("Run \#{0}, N={1}".format(n, len(stars)))
        pyplot.show()


def main(open_path, N, save_path, t_end, save, Rvir, distance, nruns, time):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    if time == 1:
        time_ev(open_path, N, save_path, t_end, save, nruns)
    else:
        #stars(open_path, N, save_path, t_end, save, nruns)
        all_runs(open_path, N, save_path, t_end, save, nruns)


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
    result.add_option("-t", dest="time", type="int", default=0,
                      help="make movie? [%default]")
    return result


if __name__ == '__main__':
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)

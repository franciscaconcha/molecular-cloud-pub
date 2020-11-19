import numpy
import matplotlib
from matplotlib import pyplot
import matplotlib.colors
from collections import defaultdict
from sklearn.neighbors import KDTree
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
import os

from amuse.lab import *
from amuse import io

from mycolors import *
from legends import *

#movie command
#"ffmpeg -framerate 5 -i {0}/%01d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {0}/movie.mp4


def local_density_vs_disk_mass(open_path, save_path, t_end, nruns, save):
    fig, (ax1, ax2) = pyplot.subplots(nrows=1,
                                            ncols=2,
                                            sharey=True,
                                            figsize=(18, 8))

    for n in range(nruns):
        path = '{0}/{1}/disks/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        # Indices for time in which star formation ends for each run
        sf_end_indices = [85, 86, 224, 127, 196, 171]

        first_file = files[sf_end_indices[n]]
        last_file = files[-1]

        min_mass = 0.1 * (0.08 | units.MSun).value_in(units.MJupiter)
        max_mass = 0.1 * (1.9 | units.MSun).value_in(units.MJupiter)

        print min_mass, max_mass

        # First file
        stars = io.read_set_from_file(path + first_file, 'hdf5', close_file=True)
        stars = stars[stars.disked]

        positions = []
        for s in stars:
            positions.append(numpy.array([s.x.value_in(units.parsec),
                                          s.y.value_in(units.parsec),
                                          s.z.value_in(units.parsec)]))
        positions = numpy.array(positions)

        tree = KDTree(positions)

        nearest_dist, nearest_ind = tree.query(positions, k=5)

        for i in range(len(stars)):
            s = stars[i]
            distances_to_neighbours = nearest_dist[i]
            max_distance = max(distances_to_neighbours)
            s.local_density = 5. / ((4./3.) * numpy.pi * max_distance**3)
            #s.local_density = 5. / (numpy.pi * max_distance ** 2)

        pyplot.set_cmap('viridis_r')

        ax1.scatter(stars.local_density,
                       stars.disk_mass.value_in(units.MJupiter),
                       s=2 * stars.disk_radius.value_in(units.au),
                       #c=stars.initial_disk_mass.value_in(units.MJupiter),
                       c='gray',
                       alpha=0.3,
                       edgecolors='none',
                       #norm=matplotlib.colors.LogNorm(),
                       #vmin=min_mass,
                       #vmax=max_mass,
                    )

        # Last file
        stars = io.read_set_from_file(path + last_file, 'hdf5', close_file=True)
        stars = stars[stars.disked]

        positions = []
        for s in stars:
            positions.append(numpy.array([s.x.value_in(units.parsec),
                                          s.y.value_in(units.parsec),
                                          s.z.value_in(units.parsec)]))
        positions = numpy.array(positions)

        tree = KDTree(positions)

        nearest_dist, nearest_ind = tree.query(positions, k=5)

        for i in range(len(stars)):
            s = stars[i]
            distances_to_neighbours = nearest_dist[i]
            max_distance = max(distances_to_neighbours)
            s.local_density = 5. / ((4. / 3.) * numpy.pi * max_distance ** 3)
            # s.local_density = 5. / (numpy.pi * max_distance ** 2)

        p = ax2.scatter(stars.local_density,
                   stars.disk_mass.value_in(units.MJupiter),
                   s=2 * stars.disk_radius.value_in(units.au),
                   #c=stars.initial_disk_mass.value_in(units.MJupiter),
                   c='gray',
                   alpha=0.3,
                   edgecolors='none',
                   #norm=matplotlib.colors.LogNorm(),
                   #vmin=min_mass,
                   #vmax=max_mass,
                    )

    #divider = make_axes_locatable(ax2)
    #cax = divider.append_axes('right', size='5%', pad=0.1)
    #cbar = fig.colorbar(p, cax=cax, orientation='vertical')

    ax1.set_xlim(left=0, right=1E7)
    ax1.set_ylim(bottom=-0.5)

    ax2.set_xlim(left=0, right=1E7)
    ax2.set_ylim(bottom=-0.5)
    ax2.yaxis.set_tick_params(labelleft=True)

    ax1.set_xscale('symlog')
    ax2.set_xscale('symlog')
    #pyplot.yscale('log')

    ax1.set_xlabel(r'Local number density (5-NN) [pc$^{-3}$]', fontsize=24)
    ax2.set_xlabel(r'Local number density (5-NN) [pc$^{-3}$]', fontsize=24)
    ax1.set_ylabel(r'Disc mass [$\mathrm{M}_{Jup}$]', fontsize=24)
    ax2.set_ylabel(r'Disc mass [$\mathrm{M}_{Jup}$]', fontsize=24)

    ax1.set_title('End of star formation')
    ax2.set_title('End of simulations')

    pyplot.tight_layout()

    if save:
        pyplot.savefig('{0}/scatterdisks.png'.format(save_path))


def distance(star1,
             star2):
    """ Return distance between star1 and star2

    :param star1: AMUSE particle
    :param star2: AMUSE particle
    :return: distance in units.parsec
    """
    return numpy.sqrt((star2.x - star1.x)**2 + (star2.y - star1.y)**2 + (star2.z - star1.z)**2)


def distance_vs_disk_mass(open_path, save_path, t_end, nruns, save):
    fig = pyplot.figure()

    for n in range(nruns):
        pyplot.clf()
        N = Ns[n]
        path = '{0}/{1}/disks/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        first_file = files[0]
        last_file = files[-1]

        stars = read_set_from_file(path + last_file, 'hdf5', close_file=True)

        #disked_stars = stars[stars.stellar_mass <= 1.9 | units.MSun]
        disked_stars = stars[stars.disked]
        bright_stars = stars[stars.bright]

        for s in disked_stars:
            distances = []
            for b in bright_stars:
                distances.append(distance(b, s).value_in(units.parsec))
            s.mean_distance = numpy.mean(distances)

        pyplot.scatter(disked_stars.mean_distance,
                       disked_stars.disk_mass.value_in(units.MJupiter),
                       s=2 * disked_stars.disk_radius.value_in(units.au),
                       c=runcolors[n],
                       alpha=0.5)#, norm=matplotlib.colors.LogNorm())

        pyplot.xscale('log')
        pyplot.yscale('log')

        pyplot.xlabel(r'Mean distance to massive stars [pc]')
        pyplot.ylabel(r'Disc mass [$\mathrm{M}_{Jup}$]')
        pyplot.title('Run \#{0} at {1:.1f} Myr'.format(n, t_end))

        if save:
            pyplot.savefig('{0}/density_vs_disk_mass/n{1}_{2:.1f}Myr.png'.format(save_path,
                                                                                 n, t_end))
        else:
            pyplot.show()


def main(open_path, N, save_path, t_end, save, Rvir, distance, nruns, movie):

    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    local_density_vs_disk_mass(open_path, save_path, t_end, nruns, save)
    #distance_vs_disk_mass(open_path, save_path, t_end, nruns, save)

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

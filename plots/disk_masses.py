import numpy
from matplotlib import pyplot
from scipy import stats
from sklearn.neighbors import KDTree

from amuse.lab import *
from amuse import io

#movie command
#"ffmpeg -framerate 5 -i {0}/%01d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {0}/movie.mp4

from legends import *  # My own custom legend definitions
from mycolors import *


def dust_mass_vs_local_density(open_path, save_path, t_end, nruns, save):
    """ Figure 6: Binned mean disc mass versus local stellar number density, projected in two dimensions.

    :param open_path: path to open results
    :param save_path: path to save figure
    :param t_end: final time to plot
    :param N: number of stars in results
    :param nruns: number of runs to plot
    :param save: if True, save figure in save_path
    """
    axs = [None, None]
    fig1, axs[0] = pyplot.subplots(1)

    #max_mass = 0.0

    #dt = 0.005
    #times = numpy.arange(0.0, t_end + dt, dt)

    for n in range(nruns):
        N = Ns[n]
        f = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, n, N, t_end)
        stars = io.read_set_from_file(f, 'hdf5', close_file=True)
        stars = stars[stars.born]
        stars = stars[stars.disked]

        all_binned_means = []
        all_binned_stds = []
        all_binned_means_locdens = []
        all_binned_stds_locdens = []

        # Calculating projected local density
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

        # Sorting disk masses by local density
        disk_masses = stars.disk_dust_mass.value_in(units.MEarth)# / 100
        masses_sorted_by_local_dens = [float(x) for _, x in sorted(zip(stars.local_density, disk_masses))]
        sorted_local_dens = sorted(stars.local_density)

        # Calculating binned masses
        binned_means_locdens = []
        binned_stds_locdens = []
        locdens_means = []
        d = 100
        for i in range(len(masses_sorted_by_local_dens)):
            if len(masses_sorted_by_local_dens) - (i + d) > 0:
                binned_means_locdens.append(numpy.median(masses_sorted_by_local_dens[i:i+d]))
                binned_stds_locdens.append(stats.sem(masses_sorted_by_local_dens[i:i+d]))
                locdens_means.append(numpy.median(sorted_local_dens[i:i+d]))
            else:
                #print "end"
                binned_means_locdens.append(numpy.median(masses_sorted_by_local_dens[i:]))
                binned_stds_locdens.append(stats.sem(masses_sorted_by_local_dens[i:]))
                locdens_means.append(numpy.median(sorted_local_dens[i:i+d]))
                break
        all_binned_means_locdens.append(numpy.array(binned_means_locdens))
        all_binned_stds_locdens.append(binned_stds_locdens)

        try:
            all_means = numpy.median(all_binned_means_locdens, axis=0)
            devs = numpy.median(all_binned_stds_locdens, axis=0)
        except:
            max_len = 0
            for a in all_binned_means_locdens:
                if len(a) > max_len:
                    max_len = len(a)

            new_sorted = []
            new_stds = []
            for a in all_binned_means_locdens:
                b = numpy.pad(a, (max_len - len(a), min(a)), 'constant')
                # constant_values=(min([min(r) for r in all_initial])))
                new_sorted.append(b)
            for a in all_binned_stds_locdens:
                b = numpy.pad(a, (max_len - len(a), min(a)), 'constant')
                # constant_values=(min([min(r) for r in all_initial])))
                new_stds.append(b)
            all_means = numpy.median(new_sorted, axis=0)
            devs = numpy.median(new_stds, axis=0)

        all_means_high = all_means + devs
        all_means_low = all_means - devs

        axs[0].plot(locdens_means,
                    all_means,
                    color=runcolors[n],
                    lw=3,
                    #label=labels[label])
                    )
        axs[0].fill_between(locdens_means,
                            all_means_high,
                            all_means_low,
                            facecolor=runcolors[n],
                            alpha=0.2)

    axs[0].set_xlabel(r'Local stellar number density [pc$^{-3}$]')
    axs[0].set_ylabel(r'Binned mean disc dust mass [$\mathrm{M}_{\oplus}$]')

    pyplot.xscale('log')
    #pyplot.yscale('log')

    if save:
        pyplot.savefig('{0}/2D_dustmass_localdensity.png'.format(save_path))
    else:
        pyplot.show()


def main(open_path, N, save_path, t_end, save, nruns):

    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    dust_mass_vs_local_density(open_path, save_path, t_end, nruns, save)

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

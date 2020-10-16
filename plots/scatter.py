import numpy
import matplotlib
from matplotlib import pyplot
import matplotlib.colors
from collections import defaultdict
from sklearn.neighbors import KDTree
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator

from amuse.lab import *
from amuse import io

from legends import *

#movie command
#"ffmpeg -framerate 5 -i {0}/%01d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {0}/movie.mp4


def local_density_vs_disk_mass_vs_initial_density(open_path, save_path, t_end, N, nruns, save, movie, log=False):
    fig = pyplot.figure()

    label = open_path.split('/')[-2].split('_')[1]
    min_dens, max_dens = 1E10, 0  # To properly set up min and max values later

    for n in range(nruns):
        # Initial stars, to calculate initial local density
        f0 = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, n, N, 0.0)
        stars0 = io.read_set_from_file(f0, 'hdf5', close_file=True)
        stars0 = stars0[stars0.disked == True]

        # Calculating initial local densities, saving as parameter of *final* stars
        positions = []
        for s in stars0:
            positions.append(numpy.array([s.x.value_in(units.parsec),
                                          s.y.value_in(units.parsec),
                                          s.z.value_in(units.parsec)]))

        positions = numpy.array(positions)
        tree = KDTree(positions)
        nearest_dist, nearest_ind = tree.query(positions, k=5)

        # Final stars
        f = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, n, N, t_end)
        stars = io.read_set_from_file(f, 'hdf5', close_file=True)
        stars = stars[stars.disked == True]

        for i in range(len(stars0)):
            s0 = stars0[i]
            s = stars[stars.key == s0.key]
            distances_to_neighbours = nearest_dist[i]
            max_distance = max(distances_to_neighbours)
            local_density = 5. / ((4./3.) * numpy.pi * max_distance**3)
            s.initial_local_density = local_density
            if local_density < min_dens:
                min_dens = local_density
            if local_density > max_dens:
                max_dens = local_density

        # Calculating final local densities
        positions = []
        for s in stars:
            positions.append(numpy.array([s.x.value_in(units.parsec),
                                          s.y.value_in(units.parsec),
                                          s.z.value_in(units.parsec)]))
        positions = numpy.array(positions)

        tree = KDTree(positions)

        nearest_dist, nearest_ind = tree.query(positions, k=5)
        #print(nearest_dist)  # drop id; assumes sorted -> see args!
        #print(nearest_ind)

        for i in range(len(stars)):
            s = stars[i]
            distances_to_neighbours = nearest_dist[i]
            max_distance = max(distances_to_neighbours)
            s.local_density = 5. / ((4./3.) * numpy.pi * max_distance**3)

        pyplot.set_cmap('viridis_r')

        pyplot.scatter(stars.local_density,
                       stars.disk_mass.value_in(units.MJupiter),
                       s=2 * stars.disk_radius.value_in(units.au),
                       c=stars.initial_local_density,
                       alpha=0.5, norm=matplotlib.colors.LogNorm())

    pyplot.xscale('symlog')
    pyplot.yscale('symlog')

    right_limit = right_limits[label]

    pyplot.gca().set_xlim(left=-0.1, right=right_limit)
    pyplot.gca().set_ylim(bottom=-0.1)

    fig.canvas.draw()

    xtick_labels = [item.get_text() for item in pyplot.gca().get_xticklabels()]
    xtick_labels[0] = r'$d_{min}$'
    ytick_labels = [item.get_text() for item in pyplot.gca().get_yticklabels()]
    ytick_labels[0] = r'$M_{min}$'

    pyplot.gca().set_xticklabels(xtick_labels)
    pyplot.gca().set_yticklabels(ytick_labels)

    cbar = pyplot.colorbar()
    pyplot.clim(min_dens, max_dens)
    cbar.set_label(r'Initial local number density [pc$^{-3}$]')

    #pyplot.legend(loc='best', fontsize=22, framealpha=1.)
    pyplot.xlabel(r'Local number density (5-NN) [pc$^{-3}$]')
    pyplot.ylabel(r'Disc mass [$\mathrm{M}_{Jup}$]')
    print label
    #print numpy.sqrt(numpy.abs(2 * stars.potential_energy().in_(units.parsec/units.Myr)))
    pyplot.suptitle(r'N = {0}, {1}, t={2:.3f} Myr'.format(N, labels[label], t_end))

    if save and movie:
        times = list(numpy.arange(0.0, t_end + 0.005, 0.005))
        print times.index(t_end)
        pyplot.savefig('{0}/{1}.png'.format(save_path, times.index(t_end)))
    elif save:
        pyplot.savefig('{0}/local_density_vs_mass_vs_ini_dens_N{1}_{2}.png'.format(save_path, N, label))


def byradius_local_density_vs_disk_mass_vs_initial_density(open_path, save_path, t_end, N, nruns, save, movie, log=False):
    times = [0.0, 1.0, 2.0]

    for folder in folders:
        fig, axs = pyplot.subplots(1, 3,
                                   figsize=(18, 6),
                                   #sharex=True,
                                   sharey=True,
                                   squeeze=True,
                                   gridspec_kw={'wspace': 0.0, 'hspace': 0.0})
        path = '{0}/{1}/'.format(open_path, folder)
        label = path.split('/')[-2].split('_')[1]

        min_dens, max_dens = 1E10, 0  # To properly set up min and max values later

        for t in times:
            for n in range(nruns):
                # Initial stars, to calculate initial local density
                f0 = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(path, n, N, 0.0)
                stars0 = io.read_set_from_file(f0, 'hdf5', close_file=True)
                stars0 = stars0[stars0.disked == True]

                # Calculating initial local densities, saving as parameter of *final* stars
                positions = []
                for s in stars0:
                    positions.append(numpy.array([s.x.value_in(units.parsec),
                                                  s.y.value_in(units.parsec),
                                                  s.z.value_in(units.parsec)]))

                positions = numpy.array(positions)
                tree = KDTree(positions)
                nearest_dist, nearest_ind = tree.query(positions, k=5)

                # Final stars
                f = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(path, n, N, t)
                stars = io.read_set_from_file(f, 'hdf5', close_file=True)
                stars = stars[stars.disked == True]

                for i in range(len(stars0)):
                    s0 = stars0[i]
                    s = stars[stars.key == s0.key]
                    distances_to_neighbours = nearest_dist[i]
                    max_distance = max(distances_to_neighbours)
                    local_density = 5. / ((4./3.) * numpy.pi * max_distance**3)
                    s.initial_local_density = local_density

                if local_density < min_dens:
                    min_dens = local_density
                if local_density > max_dens:
                    max_dens = local_density

                # Calculating final local densities
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

                pyplot.set_cmap('viridis_r')

                if t == 0.0:
                    p0 = axs[0].scatter(stars.local_density,
                                            stars.disk_mass.value_in(units.MJupiter),
                                            s=2 * stars.disk_radius.value_in(units.au),
                                            c=stars.initial_local_density,
                                            alpha=0.5,
                                            #vmin=min_dens, vmax=max_dens,
                                            norm=matplotlib.colors.LogNorm())
                elif t == 1.0:
                    p1 = axs[1].scatter(stars.local_density,
                                            stars.disk_mass.value_in(units.MJupiter),
                                            s=2 * stars.disk_radius.value_in(units.au),
                                            c=stars.initial_local_density,
                                            alpha=0.5,
                                            #vmin=min_dens, vmax=max_dens,
                                            norm=matplotlib.colors.LogNorm())
                else:  # t == 2.0
                    p2 = axs[2].scatter(stars.local_density,
                                            stars.disk_mass.value_in(units.MJupiter),
                                            s=2 * stars.disk_radius.value_in(units.au),
                                            c=stars.initial_local_density,
                                            alpha=0.5,
                                            #vmin=min_dens, vmax=max_dens,
                                            norm=matplotlib.colors.LogNorm())


        axs[0].set_xscale('symlog')
        axs[0].set_yscale('symlog')
        axs[1].set_xscale('symlog')
        axs[1].set_yscale('symlog')
        axs[2].set_xscale('symlog')
        axs[2].set_yscale('symlog')

        right_limit = right_limits[label]
        xmin = min(stars.local_density)
        ymin = min(stars.disk_mass.value_in(units.MJupiter))

        axs[0].set_xlim(left=xmin-0.2, right=right_limit)
        axs[0].set_ylim(bottom=ymin-0.2)
        axs[1].set_xlim(left=xmin-0.2, right=right_limit)
        axs[1].set_ylim(bottom=ymin-0.2)
        axs[2].set_xlim(left=xmin-0.2, right=right_limit)
        axs[2].set_ylim(bottom=ymin-0.2)

        axs[0].set_title('t = 0.0 Myr', fontsize=26)
        axs[1].set_title('t = 1.0 Myr', fontsize=26)
        axs[2].set_title('t = 2.0 Myr', fontsize=26)

        axs[0].text(0.2, 0.05,
                    labels[label],
                    ha='center', va='center', fontsize=26,
                    transform=axs[0].transAxes)

        # Fixing ticks and labels so that they don't overlap
        fig.canvas.draw()

        xtick_labels = [item.get_text() for item in axs[0].get_xticklabels()]
        xtick_labels[0] = xtick_labels[-1] = ""
        axs[0].set_xticklabels(xtick_labels)
        xticks = axs[0].xaxis.get_major_ticks()
        xticks[0].set_visible(False)
        xticks = axs[0].xaxis.get_minor_ticks()
        xticks[0].set_visible(False)

        xtick_labels = [item.get_text() for item in axs[1].get_xticklabels()]
        xtick_labels[0] = xtick_labels[-1] = ""
        axs[1].set_xticklabels(xtick_labels)
        xticks = axs[1].xaxis.get_major_ticks()
        xticks[0].set_visible(False)
        xticks = axs[1].xaxis.get_minor_ticks()
        xticks[0].set_visible(False)

        xtick_labels = [item.get_text() for item in axs[2].get_xticklabels()]
        xtick_labels[0] = ""
        axs[2].set_xticklabels(xtick_labels)
        xticks = axs[2].xaxis.get_major_ticks()
        xticks[0].set_visible(False)
        xticks = axs[2].xaxis.get_minor_ticks()
        xticks[0].set_visible(False)

        # Adding colorbar in its own new axis
        fig.subplots_adjust(top=0.9, bottom=0.2, wspace=0.33)

        divider = make_axes_locatable(axs[2])
        cax = divider.append_axes('right', size='5%', pad=0.1)
        cbar = fig.colorbar(p2, cax=cax, orientation='vertical')

        cbar.set_label(r'Initial local number density [pc$^{-3}$]')

        # x-label
        fig.text(0.5, 0.05,
                 r'Local stellar number density [pc$^{-3}$]',
                 ha='center', va='center', fontsize=26)

        # y-label
        fig.text(0.08, 0.55,
                 r'Total disc mass [$\mathrm{M}_{Jup}$]',
                 ha='center', va='center', rotation='vertical', fontsize=26)

        if save:
            pyplot.savefig('{0}/discmass_vs_localdens_{1}.png'.format(save_path, label))


def local_density_vs_disk_mass_vs_stellar_velocity(open_path, save_path, t_end, N, nruns, save, movie,
                                                   min_cb=0.1, max_cb=50., log=False):
    fig = pyplot.figure()

    label = open_path.split('/')[-2].split('_')[1]

    # To assign colorbar limits
    min_vel, min_escape_vel = 1E10, 1E10  # To properly set up min and max values later
    all_densities = []
    all_masses = []
    all_radii = []

    for n in range(nruns):
        f = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, n, N, t_end)
        all_stars = io.read_set_from_file(f, 'hdf5', close_file=True)
        stars = all_stars[all_stars.disked == True]

        positions = []
        for s in stars:
            positions.append(numpy.array([s.x.value_in(units.parsec),
                                          s.y.value_in(units.parsec),
                                          s.z.value_in(units.parsec)]))
        positions = numpy.array(positions)

        tree = KDTree(positions)

        nearest_dist, nearest_ind = tree.query(positions, k=5)
        #print(nearest_dist)  # drop id; assumes sorted -> see args!
        #print(nearest_ind)

        for i in range(len(stars)):
            s = stars[i]
            distances_to_neighbours = nearest_dist[i]
            max_distance = max(distances_to_neighbours)
            s.local_density = 5. / ((4./3.) * numpy.pi * max_distance**3)
        #print local_densities

            s.velocity_mag = numpy.sqrt(s.vx.value_in(units.km/units.s)**2
                                        + s.vy.value_in(units.km/units.s)**2
                                        + s.vz.value_in(units.km/units.s)**2)

        escape_vel2 = (2 * constants.G * all_stars.stellar_mass.sum()) / all_stars.virial_radius()
        escape_vel = numpy.sqrt(escape_vel2.value_in(units.km**2/units.s**2))

        print escape_vel
        print len(stars.velocity_mag[stars.velocity_mag >= escape_vel])

        if escape_vel <= min_escape_vel:
            min_escape_vel = escape_vel
        if min(stars.velocity_mag) <= min_vel:
            min_vel = min(stars.velocity_mag)

        pyplot.set_cmap('YlGnBu')

        all_densities.append(stars.local_density)
        all_masses.append(stars.disk_mass.value_in(units.MJupiter))
        all_radii.append(2 * stars.disk_radius.value_in(units.au))

        p = pyplot.scatter(stars.local_density,
                       stars.disk_mass.value_in(units.MJupiter),
                       s=2 * stars.disk_radius.value_in(units.au),
                       c=stars.velocity_mag,
                       alpha=0.5, norm=matplotlib.colors.LogNorm())

        pyplot.gca().scatter(stars[stars.velocity_mag > escape_vel].local_density,
                       stars[stars.velocity_mag > escape_vel].disk_mass.value_in(units.MJupiter),
                       s=2 * stars[stars.velocity_mag > escape_vel].disk_radius.value_in(units.au),
                       c='r',
                       alpha=0.5)#, norm=matplotlib.colors.LogNorm())

    pyplot.xscale('symlog')
    pyplot.yscale('symlog')

    right_limit = right_limits[label]

    pyplot.gca().set_xlim(left=-0.1, right=right_limit)
    pyplot.gca().set_ylim(bottom=-0.1, top=250)

    fig.canvas.draw()

    xtick_labels = [item.get_text() for item in pyplot.gca().get_xticklabels()]
    xtick_labels[0] = r'$d_{min}$'
    ytick_labels = [item.get_text() for item in pyplot.gca().get_yticklabels()]
    ytick_labels[0] = r'$M_{min}$'

    pyplot.gca().set_xticklabels(xtick_labels)
    pyplot.gca().set_yticklabels(ytick_labels)

    cbar = pyplot.colorbar(p)

    cbar.set_label(r'Stellar velocity [km/s]')

    if nruns > 0:
        pyplot.clim([min_cb, max_cb])
    else:
        pyplot.clim([min_vel, min_escape_vel])

    #pyplot.legend(loc='best', fontsize=22, framealpha=1.)
    pyplot.xlabel(r'Local number density (5-NN) [pc$^{-3}$]')
    pyplot.ylabel(r'Disc mass [$\mathrm{M}_{Jup}$]')
    print label
    #print numpy.sqrt(numpy.abs(2 * stars.potential_energy().in_(units.parsec/units.Myr)))
    pyplot.suptitle(r'N = {0}, {1}, t={2:.3f} Myr'.format(N, labels[label], t_end))

    if save and movie:
        times = list(numpy.arange(0.0, t_end + 0.005, 0.005))
        print times.index(t_end)
        pyplot.savefig('{0}/{1}.png'.format(save_path, times.index(t_end)))
    elif save:
        pyplot.savefig('{0}/local_density_vs_mass_vs_stellar_vel_N{1}_{2}.png'.format(save_path, N, label))

def main(open_path, N, save_path, t_end, save, Rvir, distance, nruns, movie):

    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    paths = ['results/large_cluster/N1E3_R01',
             'results/large_cluster/N1E3_R03',
             'results/large_cluster/N1E3_R05',
             'results/large_cluster/N1E3_R1',
             'results/large_cluster/N1E3_R25',
             'results/large_cluster/N1E3_R5']

    #local_density_vs_disk_mass_vs_initial_density(open_path, save_path, t_end, N, nruns, save, movie, log=False)
    #local_density_vs_disk_mass_vs_stellar_velocity(open_path, save_path, t_end, N, nruns, save, movie, log=False)
    #local_density_vs_disk_mass_tracks(open_path, save_path, t_end, N, nruns, save, log=False)

    byradius_local_density_vs_disk_mass_vs_initial_density(open_path, save_path, t_end, N, nruns, save, movie, log=False)

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

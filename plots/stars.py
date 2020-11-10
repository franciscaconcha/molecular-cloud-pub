import numpy
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

from amuse.lab import *
from amuse import io
from amuse.ext.basicgraph import Graph, MinimumSpanningTreeFromEdges, MinimumSpanningTree

from mycolors import *
from legends import *

# movie command
# "ffmpeg -framerate 5 -i {0}/%01d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {0}/movie.mp4


def stars(open_path, i, save_path, save, nrun):
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

    #f = '{0}/{1}/gravity_stars.hdf5'.format(open_path, nrun)
    path = '{0}/{1}/prep/'.format(open_path, nrun)
    files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
    files = [x for x in files if '.hdf5' in x]
    files.sort(key=lambda f: float(filter(str.isdigit, f)))
    last_file = files[i]

    print last_file

    stars = read_set_from_file(path + last_file, 'hdf5', close_file=True)
    t = float(last_file.split('t')[1].split('.hdf5')[0])

    disked_stars = stars[stars.disked]
    massive_stars = stars[stars.bright]

    if t == 1.724:
        print t
        for s in disked_stars:
            print s.disk_radius.in_(units.au), s.initial_disk_radius.in_(units.au), s.disk_mass.in_(units.MJupiter), s.stellar_mass.in_(units.MSun), s.tborn.in_(units.Myr)
    elif t == 1.766:
        print t
        for s in disked_stars:
            print s.disk_radius.in_(units.au), s.initial_disk_radius.in_(units.au), s.disk_mass.in_(units.MJupiter), s.stellar_mass.in_(units.MSun), s.tborn.in_(units.Myr)
    elif t == 1.808:
        print t
        for s in disked_stars:
            print s.disk_radius.in_(units.au), s.initial_disk_radius.in_(units.au), s.disk_mass.in_(units.MJupiter), s.stellar_mass.in_(units.MSun), s.tborn.in_(units.Myr)
        return


    pyplot.set_cmap('viridis_r')

    s = ax.scatter(disked_stars.x.value_in(units.parsec),
               disked_stars.y.value_in(units.parsec),
               marker='o',
               s=disked_stars.disk_radius.value_in(units.au),
               c=disked_stars.disk_mass.value_in(units.MJupiter),
               alpha=0.5,
               #norm=matplotlib.colors.LogNorm())
                   )

    ax.scatter(massive_stars.x.value_in(units.parsec),
               massive_stars.y.value_in(units.parsec),
               marker='*',
               s=100,
               c='red',
               alpha=0.5,
               lw=1)

    ax.set_xlabel("x [pc]")
    ax.set_ylabel("y [pc]")

    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])

    ax.set_aspect('equal')

    #divider = make_axes_locatable(axs[2])
    #cax = divider.append_axes('right', size='5%', pad=0.1)
    cbar = fig.colorbar(s)#p2, cax=cax, orientation='vertical')

    cbar.set_label(r'Initial local number density [pc$^{-3}$]')

    fig.suptitle("Run \#{0}, N={1}, t = {2} Myr".format(nrun,
                                                        len(stars),
                                                        t))

    if save:
        pyplot.savefig('{0}/{1}.png'.format(save_path, i))
    else:
        pyplot.show()


def mst(open_path, N, save_path, t_end, save, nrun):
    fig = pyplot.figure()
    ax = pyplot.gca()

    f = '{0}/{1}/gravity_stars.hdf5'.format(open_path, nrun)
    stars = io.read_set_from_file(f, 'hdf5', close_file=True)
    stars.mass = stars.stellar_mass

    N = len(stars)

    graph = Graph()

    def distfunc(p, q):
        return (((p.x - q.x) ** 2 + (p.y - q.y) ** 2) ** 0.5).value_in(p.x.unit)

    from amuse.units.quantities import zero

    for p in stars:
        d = distfunc(p, stars)
        for i, q in enumerate(stars):
            if p != q:
                graph.add_edge(p, q, d[i])

    all_edges = graph.all_edges()

    ml = reduce(lambda x, y: x + y[0], all_edges, zero) / len(all_edges)

    mst = MinimumSpanningTreeFromEdges(all_edges)

    mlmst = reduce(lambda x, y: x + y[0], mst, zero) / len(mst)
    # normalize
    mlmst = mlmst / (N * numpy.pi) ** 0.5 * (N - 1)
    print mlmst

def all_runs(open_path, N, save_path, t_end, save, nruns):
    fig, axs = pyplot.subplots(2, 5,
                               figsize=(16, 8),
                               subplot_kw=dict(aspect='equal',
                                               adjustable='box-forced'))  # rows, columns
    fig.subplots_adjust(wspace=0.5,
                        hspace=0.05)

    ax = {0: axs[0, 0],
          1: axs[0, 1],
          2: axs[0, 2],
          3: axs[0, 3],
          4: axs[0, 4],
          5: axs[1, 0],
          6: axs[1, 1],
          7: axs[1, 2],
          8: axs[1, 3],
          9: axs[1, 4],
          #10: axs[3, 1],
          #11: axs[3, 2]
    }

    for n in range(nruns):
        path = '{0}/{1}/'.format(open_path, n)
        #files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        #star_files = [x for x in files if '.hdf5' in x]
        #star_files.sort(key=lambda f: int(filter(str.isdigit, f)))
        f = '{0}/{1}/gravity_stars.hdf5'.format(open_path, n)
        stars = io.read_set_from_file(f, 'hdf5', close_file=True)

        #print stars[0].x.in_(units.parsec), stars[0].y.in_(units.parsec), stars[0].z.in_(units.parsec)

        #f = '{0}/{1}/disks/0/{2}'.format(open_path, n, star_files[0])
        #stars = io.read_set_from_file(f, 'hdf5', close_file=True)

        #print stars[0].x.in_(units.parsec), stars[0].y.in_(units.parsec), stars[0].z.in_(units.parsec)

        disked_stars = stars[stars.stellar_mass <= 1.9 | units.MSun]
        massive_stars = stars[stars.stellar_mass > 1.9 | units.MSun]

        #print n + 1, len(stars), len(stars[stars.born])

        ax[n].scatter(disked_stars.x.value_in(units.parsec),
                       disked_stars.y.value_in(units.parsec),
                       marker='o',
                       s=stars.disk_radius.value_in(units.au),
                       c=runcolors[n],
                       alpha=0.5,
                       lw=1)

        ax[n].scatter(massive_stars.x.value_in(units.parsec),
                       massive_stars.y.value_in(units.parsec),
                       marker='*',
                       #s=stars.disk_radius.value_in(units.au),
                       c='k',
                       alpha=0.5,
                       lw=1)

        #ax0.set_xlabel("x [pc]")
        #ax0.set_ylabel("y [pc]")

        ax[n].set_xlim([-2, 2])
        ax[n].set_ylim([-2, 2])

    ax[2].set_xlabel("x [pc]")
    ax[7].set_xlabel("x [pc]")

    ax[0].set_ylabel("y [pc]")
    ax[5].set_ylabel("y [pc]")

    #pyplot.tight_layout()

    if save:
        pyplot.savefig('{0}/stars.png'.format(save_path))
    else:
        pyplot.show()



def main(open_path, i, N, save_path, t_end, save, movie, nrun):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    if movie:
        path = '{0}/{1}/prep/'.format(open_path, nrun)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        l = len(files)
        for i in range(l):
            stars(open_path, i, save_path, 1, nrun)
            pyplot.close('all')
    else:
        stars(open_path, i, save_path, save, nrun)

    #mst(open_path, N, save_path, t_end, save, nruns)

    #all_runs(open_path, N, save_path, t_end, save, nruns)


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()

    result.add_option("-p", dest="open_path", type="string", default='/media/fran/data1/photoevap/results',
                      help="path to results to plot [%default]")
    result.add_option("-i", dest="i", type="int", default=1,
                      help="index of file to plot [%default]")
    result.add_option("-N", dest="N", type="int", default=1000,
                      help="number of stars [%default]")
    result.add_option("-S", dest="save", type="int", default=0,
                      help="save plot? [%default]")
    result.add_option("-s", dest="save_path", type="string", default='/media/fran/data1/photoevap-paper/figures',
                      help="path to save the results [%default]")
    result.add_option("-e", dest="t_end", type="float", default='5.0',
                      help="end time to use for plots [%default]")
    result.add_option("-m", dest="movie", type="int", default=0,
                      help="make movie? [%default]")
    result.add_option("-n", dest="nrun", type="int", default=1,
                      help="number of runs to plot for averages [%default]")
    return result


if __name__ == '__main__':
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)

import numpy
from matplotlib import pyplot
from scipy import stats
from sklearn.neighbors import KDTree

from amuse.lab import *
from amuse import io

#movie command
#"ffmpeg -framerate 5 -i {0}/%01d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {0}/movie.mp4

from legends import *  # My own custom legend definitions


def projected_mass_vs_local_density(open_path, save_path, t_end, N, nruns, save):
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

    max_mass = 0.0

    dt = 0.005
    times = numpy.arange(0.0, t_end + dt, dt)

    for n in range(nruns):
        N = Ns[n]
        f = '{0}/{1}/N{2}_t{3:.3f}.hdf5'.format(open_path, n, N, t_end)
        stars = io.read_set_from_file(f, 'hdf5', close_file=True)

        all_binned_means = []
        all_binned_stds = []
        all_binned_means_locdens = []
        all_binned_stds_locdens = []

        positions = []
        for s in stars:
            positions.append(numpy.array([s.x.value_in(units.parsec),
                                          s.y.value_in(units.parsec)]))#,
                                          #s.z.value_in(units.parsec)]))
        positions = numpy.array(positions)

        tree = KDTree(positions)

        nearest_dist, nearest_ind = tree.query(positions, k=5)
        # print(nearest_dist)  # drop id; assumes sorted -> see args!
        # print(nearest_ind)

        for i in range(len(stars)):
            s = stars[i]
            distances_to_neighbours = nearest_dist[i]
            max_distance = max(distances_to_neighbours)
            s.local_density = 5. / (numpy.pi * max_distance ** 2)

            disk_masses = stars.disk_mass.value_in(units.MEarth) / 100
            masses_sorted_by_local_dens = [float(x) for _, x in sorted(zip(stars.local_density, disk_masses))]
            sorted_local_dens = sorted(stars.local_density)

            # LOCAL DENSITY
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
                    #color=colors[label],
                    lw=3,
                    #label=labels[label])
                    )
        axs[0].fill_between(locdens_means,
                            all_means_high,
                            all_means_low,
                            #facecolor=colors[label],
                            alpha=0.2)

    """from astropy.table import Table
    data = Table.read('data/surfacedensities_var.dat', format='ascii.ecsv')
    point_ONC = [float(data[0]['surfdens_YSO']),
                 float(data[0]['Mavg']),
                 numpy.log(float(data[0]['sigma_ln_Mavg/Mearth'])),
                 data[0]['Region']]
    point_Lupus = [float(data[1]['surfdens_YSO']),
                   float(data[1]['Mavg']),
                   numpy.log(float(data[1]['sigma_ln_Mavg/Mearth'])),
                   data[1]['Region']]
    point_OMC2 = [float(data[2]['surfdens_YSO']),
                  float(data[2]['Mavg']),
                  numpy.log(float(data[2]['sigma_ln_Mavg/Mearth'])),
                  data[2]['Region']]
    point_NGCE = [float(data[3]['surfdens_YSO']),
                  float(data[3]['Mavg']),
                  numpy.log(float(data[3]['sigma_ln_Mavg/Mearth'])),
                  data[3]['Region']]
    point_NGCW = [float(data[4]['surfdens_YSO']),
                  float(data[4]['Mavg']),
                  numpy.log(float(data[4]['sigma_ln_Mavg/Mearth'])),
                  data[4]['Region']]
    point_Taurus = [float(data[5]['surfdens_YSO']),
                    float(data[5]['Mavg']),
                    numpy.log(float(data[5]['sigma_ln_Mavg/Mearth'])),
                    data[5]['Region']]

    surfdens = [point_ONC[0], point_Lupus[0], point_OMC2[0], point_NGCE[0], point_NGCW[0], point_Taurus[0]]
    obsmasses = [point_ONC[1], point_Lupus[1], point_OMC2[1], point_NGCE[1], point_NGCW[1], point_Taurus[1]]
    sigmas = [point_ONC[2], point_Lupus[2], point_OMC2[2], point_NGCE[2], point_NGCW[2], point_Taurus[2]]
    names = [point_ONC[3], point_Lupus[3], point_OMC2[3], point_NGCE[3], point_NGCW[3], point_Taurus[3]]

    #pyplot.scatter(surfdens, obsmasses, c='k', marker='x', s=80, lw=4)
    markers, caps, bars = pyplot.errorbar(surfdens[:3], obsmasses[:3], yerr=sigmas[:3],
                                          color='navy',
                                          marker='D',
                                          markersize=10,
                                          #markeredgewidth=6,
                                          elinewidth=2,
                                          capsize=2,
                                          ls="None",
                                          #alpha=0.8,
                                          zorder=10)
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(surfdens[3], obsmasses[3], yerr=sigmas[3],
                                          color='royalblue',
                                          marker='D',
                                          markersize=10,
                                          #markeredgewidth=6,
                                          elinewidth=2,
                                          capsize=2,
                                          ls="None",
                                          #alpha=0.5,
                                          zorder=12)
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(surfdens[4:], obsmasses[4:], yerr=sigmas[4:],
                                          color='navy',
                                          marker='D',
                                          markersize=10,
                                          #markeredgewidth=6,
                                          elinewidth=2,
                                          capsize=2,
                                          ls="None",
                                          #alpha=0.5,
                                          zorder=15)
    [bar.set_alpha(0.5) for bar in bars]

    # Positions for text labels
    xlocs = {data[0]['Region']: surfdens[0] * 1.25,  # ONC
             data[1]['Region']: surfdens[1] * 0.26,  # Lupus
             data[2]['Region']: surfdens[2] * 1.25,  # OMC-2
             data[3]['Region']: surfdens[3] * 1.2,  # NGC2024 East
             data[4]['Region']: surfdens[4] * 1.2,  # NGC2024 West
             data[5]['Region']: surfdens[5] * 0.22}  # Taurus

    ylocs = {data[0]['Region']: obsmasses[0] * 0.9,  # ONC
             data[1]['Region']: obsmasses[1] * 0.9,  # Lupus
             data[2]['Region']: obsmasses[2] * 0.9,  # OMC-2
             data[3]['Region']: obsmasses[3] * 0.9,  # NGC2024 East
             data[4]['Region']: obsmasses[4] * 0.9,  # NGC2024 West
             data[5]['Region']: obsmasses[5] * 0.9}  # Taurus

    for i, txt in enumerate(names):
        if txt == data[3]['Region']:
            axs[0].annotate(txt,
                            (xlocs[data[i]['Region']], ylocs[data[i]['Region']]),
                            fontsize=20,
                            color='royalblue')
        else:
            axs[0].annotate(txt,
                            (xlocs[data[i]['Region']], ylocs[data[i]['Region']]),
                            fontsize=20,
                            color='navy')"""

    axs[0].set_xlabel(r'Projected local stellar number density [pc$^{-2}$]')
    axs[0].set_ylabel(r'Binned mean disc dust mass [$\mathrm{M}_{\oplus}$]')

    first_legend = pyplot.legend([R01ShadedObject(), R03ShadedObject(), R05ShadedObject(),
                                  R1ShadedObject(), R25ShadedObject(), R5ShadedObject()],
                                  [labels['R01'], labels['R03'], labels['R05'],
                                   labels['R1'], labels['R25'], labels['R5'],
                                   ],
                                  handler_map={R01ShadedObject: R01ShadedObjectHandler(),
                                               R03ShadedObject: R03ShadedObjectHandler(),
                                               R05ShadedObject: R05ShadedObjectHandler(),
                                               R1ShadedObject: R1ShadedObjectHandler(),
                                               R25ShadedObject: R25ShadedObjectHandler(),
                                               R5ShadedObject: R5ShadedObjectHandler()},
                                  loc='lower left',
                                  ncol=2,
                                  fontsize=18, framealpha=0.4)

    pyplot.gca().add_artist(first_legend)

    pyplot.legend([DottedShadedObject(), SolidShadedObject()],
                  [r"t = 0.0 Myr",
                   r"t = 2.0 Myr",
                   ],
                  handler_map={DottedShadedObject: DottedShadedObjectHandler(),
                               SolidShadedObject: SolidShadedObjectHandler()},
                  loc='lower left',
                  bbox_to_anchor=(0.0, 0.2),
                  fontsize=18, framealpha=0.4)

    pyplot.xscale('log')
    #pyplot.yscale('log')

    if save:
        pyplot.savefig('{0}/2D_dustmass_localdensity.png'.format(save_path))


def main(open_path, N, save_path, t_end, save, nruns):

    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    projected_mass_vs_local_density(open_path, save_path, t_end, N, nruns, save)

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

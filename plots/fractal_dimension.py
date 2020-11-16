from matplotlib import pyplot
import os

from amuse.lab import *
from amuse import io

from mycolors import *


def fractal_dimension(open_path, nruns, save, save_path):
    for n in range(nruns):
        path = '{0}/{1}/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if 'hydro_stars' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        stars = read_set_from_file(path + files[-1], 'hdf5', close_file=True)
        #t = float(files[-1].split('t')[1].split('.hdf5')[0])

        # Have to do this because of some silly book keeping choices
        # I add 10 Myr to the last step of the molecular cloud collapse
        # to break out of the loop
        t = stars.get_timestamp()
        t -= 10 | units.Myr
        if t < 1.0 | units.Myr:
            t += 10 | units.Myr
        fd = stars.box_counting_dimension()

        print "Run {0}, t = {1:.3f} Myr, Fd = {2}".format(n,
                                                          t.value_in(units.Myr),
                                                          fd)


def fd_vs_time(open_path, nruns, save, save_path):
    fig1, axs1 = pyplot.subplots(1)

    times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    dimension = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

    for n in range(nruns):
        path = '{0}/{1}/disks/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        for f in files[1:]:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            t = float(f.split('t')[1].split('.hdf5')[0])
            fd = stars.box_counting_dimension()
            times[n].append(t)
            dimension[n].append(fd)

    for n in range(nruns):
        axs1.plot(times[n],
                  dimension[n],
                  lw=3,
                  c=runcolors[n],
                  label=r'Run \#{0}'.format(n))

    axs1.xlabel('Time [Myr]')
    axs1.ylabel(r'$F_d$')

    axs1.legend()

    if save:
        pyplot.savefig('{0}/Fd.png'.format(save_path))


def fd_scatter(open_path, nruns, save, save_path):

    # fractal dimensions of each run at tendSF as measured in fractal_dimension
    Fd = {0: 1.6,
          1: 1.5,
          2: 1.7,
          3: 1.5,
          4: 1.4,
          5: 1.6}

    # times of end of star formation (Myr)
    times = {0: 4.25,
             1: 4.29,
             2: 10.01,
             3: 5.93,
             4: 8.87,
             5: 7.74}

    fig = pyplot.figure()

    pyplot.scatter(x=times.values(),
                   y=Fd.values(),
                   c='k',
                   marker='D',
                   s=120,
                   alpha=0.8,
                   label='Simulations'
                   )

    # Observational points
    refs = {'CW2004': r'Cartwright \& Whitworth 2004',
            'H2002': r'Hartmann 2002',
            'KH2008': r'Kraus \& Hillenbrand 2008',
            'Simon1997': r'Simon 1997'}

    from astropy.table import Table
    data = Table.read('data/Fd_data.txt', format='ascii.ecsv')
    CW2004 = data[data['Fd_source'] == 'CW2004']
    H2002 = data[data['Fd_source'] == 'H2002']
    KH2008 = data[data['Fd_source'] == 'KH2008']
    Simon1997 = data[data['Fd_source'] == 'Simon1997']

    markers, caps, bars = pyplot.errorbar(CW2004['Age'],
                                            CW2004['Fd'],
                                            xerr=CW2004['Age_error'],
                                            yerr=CW2004['Fd_error'],
                                            c=colors['yellow'],
                                            label=refs['CW2004'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None')
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(H2002['Age'],
                                            H2002['Fd'],
                                            xerr=H2002['Age_error'],
                                            yerr=H2002['Fd_error'],
                                            c=colors['red'],
                                            label=refs['H2002'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None'
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(KH2008['Age'],
                                            KH2008['Fd'],
                                            xerr=KH2008['Age_error'],
                                            yerr=KH2008['Fd_error'],
                                            c=colors['turquoise'],
                                            label=refs['KH2008'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None'
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(Simon1997['Age'],
                                            Simon1997['Fd'],
                                            xerr=Simon1997['Age_error'],
                                            yerr=Simon1997['Fd_error'],
                                            c=colors['orange'],
                                            label=refs['Simon1997'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None'
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    #pyplot.xscale('log')
    pyplot.legend(fontsize=18)
    pyplot.xlabel('Age [Myr]')
    pyplot.ylabel(r'$F_d$', fontsize=24)


def main(open_path, nruns, save, save_path):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    #fractal_dimension(open_path, nruns, save, save_path)
    #fd_vs_time(open_path, nruns, save, save_path)
    fd_scatter(open_path, nruns, save, save_path)

    if not save:
        pyplot.show()


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()

    result.add_option("-p", dest="open_path", type="string", default='/media/fran/data1/photoevap/results',
                      help="path to results to plot [%default]")
    result.add_option("-n", dest="nruns", type="int", default=1,
                      help="number of runs to plot for averages [%default]")
    result.add_option("-S", dest="save", type="int", default=0,
                      help="save plot? [%default]")
    result.add_option("-s", dest="save_path", type="string", default='/media/fran/data1/photoevap-paper/figures',
                      help="path to save the results [%default]")
    return result


if __name__ == '__main__':
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)

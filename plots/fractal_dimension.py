from matplotlib import pyplot
import os
import numpy

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

    """for n in range(nruns):
        path = '{0}/{1}/disks/'.format(open_path, n)
        print path
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        for f in files[1:]:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            t = float(f.split('t')[1].split('.hdf5')[0])
            fd = stars.box_counting_dimension()
            times[n].append(t)
            dimension[n].append(fd)

    print times
    print dimension"""

    from read_fd import times, dimension, end_times

    # Find the index in time[n] when star formation ends
    indexes = []
    for n in range(nruns):
        i = 0
        for t in times[n]:
            if t < end_times[n]:
                i += 1
        indexes.append(i)

    for n in range(nruns):
        i = indexes[n]
        axs1.plot(times[n][:i],
                  dimension[n][:i],
                  lw=3,
                  c=runcolors[n],
                  label=r'Run \#{0}'.format(n))
        axs1.plot(times[n][i:],
                  dimension[n][i:],
                  lw=3,
                  ls=":",
                  c=runcolors[n],
                  )

    path = '{0}/plummer6k/0/'.format(open_path)
    print path
    files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
    files = [x for x in files if '.hdf5' in x]
    files.sort(key=lambda f: float(filter(str.isdigit, f)))

    plummert, plummerfd = [], []

    for f in files[1:]:
        stars = read_set_from_file(path + f, 'hdf5', close_file=True)
        t = float(f.split('t')[1].split('.hdf5')[0])
        fd = stars.box_counting_dimension()
        plummert.append(t)
        plummerfd.append(fd)

    print plummert
    print plummerfd

    axs1.plot(plummert,
              plummerfd,
              lw=3,
              c='k',
              label=r'Plummer sphere')

    axs1.set_xlabel('Time [Myr]')
    axs1.set_ylabel(r'$F_d$')

    pyplot.legend(loc='upper right', ncol=2, fontsize=20)

    if save:
        pyplot.savefig('{0}/Fd_vs_time.png'.format(save_path))


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
            'Simon1997': r'Simon 1997',
            'Parker2014': r'Parker 2014'}

    from astropy.table import Table
    data = Table.read('data/Fd_data.txt', format='ascii.ecsv')
    CW2004 = data[data['Fd_source'] == 'CW2004']
    H2002 = data[data['Fd_source'] == 'H2002']
    KH2008 = data[data['Fd_source'] == 'KH2008']
    Simon1997 = data[data['Fd_source'] == 'Simon1997']
    #Sanchez2009 = data[data['Fd_source'] == 'Sanchez2009']

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
                                            ls='None',
                                            fillstyle='none',
                                            mew=2,  # marker edge width
                                            )
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
                                            ls='None',
                                            zorder=1)
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

    """markers, caps, bars = pyplot.errorbar(numpy.exp(Sanchez2009['Age']),
                                            Sanchez2009['Fd'],
                                            xerr=Sanchez2009['Age_error'],
                                            yerr=Sanchez2009['Fd_error'],
                                            c=colors['brown'],
                                            label=refs['Sanchez2009'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None',
                                            #alpha=0.5,
                                            )
    [bar.set_alpha(0.5) for bar in bars]"""

    for i, txt in enumerate(CW2004['Region']):
        if txt == 'Taurus':
            pyplot.annotate(txt,
                                (CW2004[i]['Age'] - 1.1, CW2004[i]['Fd'] - 0.10),
                                fontsize=20,
                                color=colors['yellow'])
        else:
            pyplot.annotate(txt,
                                (CW2004[i]['Age'] - 1.8, CW2004[i]['Fd'] - 0.10),
                                fontsize=20,
                                color=colors['yellow'])

    for i, txt in enumerate(H2002['Region']):
        pyplot.annotate(txt,
                            (H2002[i]['Age'] + .1, H2002[i]['Fd'] - 0.10),
                            fontsize=20,
                            color=colors['red'])

    for i, txt in enumerate(KH2008['Region']):
        pyplot.annotate(txt,
                            (KH2008[i]['Age'] + .1, KH2008[i]['Fd'] + 0.02),
                            fontsize=20,
                            color=colors['turquoise'])

    for i, txt in enumerate(Simon1997['Region']):
        if txt == 'Ophiuchus':
            pyplot.annotate(txt,
                                (Simon1997[i]['Age'] + .1, Simon1997[i]['Fd'] + 0.02),
                                fontsize=20,
                                color=colors['orange'])
        elif txt == 'Trapezium':
            pyplot.annotate(txt,
                                (Simon1997[i]['Age'] - 1.3, Simon1997[i]['Fd'] + 0.02),
                                fontsize=16,
                                color=colors['orange'])

    #pyplot.xscale('log')
    pyplot.legend(fontsize=18)
    pyplot.xlabel('Age [Myr]')
    pyplot.ylabel(r'$F_d$', fontsize=24)

    if save:
        pyplot.savefig('{0}/Fd_vs_age.png'.format(save_path))


def main(open_path, nruns, save, save_path):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    #fractal_dimension(open_path, nruns, save, save_path)
    fd_vs_time(open_path, nruns, save, save_path)
    #fd_scatter(open_path, nruns, save, save_path)

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

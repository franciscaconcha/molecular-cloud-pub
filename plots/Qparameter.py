from matplotlib import pyplot
import os

from amuse.lab import *
from amuse import io

from mycolors import *


def Q_vs_time(open_path, nruns, save, save_path):
    fig = pyplot.figure()

    times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    Qparam = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

    for n in range(nruns):
        path = '{0}/{1}/disks/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        print path 

        for f in files[1:]:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            t = float(f.split('t')[1].split('.hdf5')[0])
            #Rvir = stars.virial_radius()
            #converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), Rvir)
            #stars = stars.bound_subset(tidal_radius=Rvir, unit_converter=converter)
            q = stars.Qparameter()
            times[n].append(t)
            Qparam[n].append(q)

    print times
    print Qparam

    for n in range(nruns):
        pyplot.plot(times[n],
                    Qparam[n],
                    lw=3,
                    c=runcolors[n],
                    label=r'Run \#{0}'.format(n))

    plummert, plummerq = [], []

    for f in files[1:]:
        stars = read_set_from_file(path + f, 'hdf5', close_file=True)
        t = float(f.split('t')[1].split('.hdf5')[0])
        #Rvir = stars.virial_radius()
        #converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), Rvir)
        #stars = stars.bound_subset(tidal_radius=Rvir, unit_converter=converter)
        q = stars.Qparameter()
        plummert.append(t)
        plummerq.append(q)

    print plummert
    print plummerq

    pyplot.plot(plummert,
              plummerq,
              lw=3,
              c='k',
              label=r'Plummer sphere')

    pyplot.xlabel('Time [Myr]')
    pyplot.ylabel('Q parameter')

    pyplot.legend()

    if save:
        pyplot.savefig('{0}/Q_vs_time.png'.format(save_path))


def Q_parameter(open_path, nruns):
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
        Rvir = stars.virial_radius()
        converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), Rvir)
        stars = stars.bound_subset(tidal_radius=Rvir, unit_converter=converter)
        Q = stars.Qparameter()

        print "Run {0}, t = {1:.3f} Myr, Q = {2}".format(n,
                                                          t.value_in(units.Myr),
                                                          Q)


def Q_scatter(open_path, nruns, save, save_path):
    # Q of each run at tendSF as measured in fractal_dimension
    Q = {0: 0.816,
         1: 0.776,
         2: 0.705,
         3: 0.745,
         4: 0.743,
         5: 0.663}

    # times of end of star formation (Myr)
    times = {0: 4.25,
             1: 4.29,
             2: 10.01,
             3: 5.93,
             4: 8.87,
             5: 7.74}

    fig = pyplot.figure()

    pyplot.scatter(x=times.values(),
                   y=Q.values(),
                   c='k',
                   marker='D',
                   s=120,
                   alpha=0.8,
                   label='Simulations'
                   )

    # Observational points
    refs = {'CW2004': r'Cartwright \& Whitworth 2004',
            'H2002': r'Hartmann 2002',
            'H1998': r'Hillenbrand \& Hartmann 1998',
            'KH2008': r'Kraus \& Hillenbrand 2008',
            'Simon1997': r'Simon 1997',
            'Parker2017': r'Parker \& Alves de Oliveira 2017',
            'Wright2013': r'Wright \emph{et al.} 2014',
            'Parker2014': 'Parker 2014',
            'N2008': r'Neuh\"auser \& Forbrich 2008'}

    from astropy.table import Table
    data = Table.read('data/Fd_data.txt', format='ascii.ecsv')
    H1998 = data[data['Q_source'] == 'Hillenbrand1998']
    CW2004 = data[data['Q_source'] == 'CW2004']
    N2008 = data[data['Q_source'] == 'Neuhauser2008']
    KH2008 = data[data['Q_source'] == 'KH2008']
    Wright2013 = data[data['Q_source'] == 'Wright2014']
    Parker2014 = data[data['Q_source'] == 'Parker2014']
    Parker2017 = data[data['Q_source'] == 'Parker2017']

    markers, caps, bars = pyplot.errorbar(H1998['Age'],
                                            H1998['Q'],
                                            xerr=H1998['Age_error'],
                                            yerr=H1998['Q_error'],
                                            c=colors['pink'],
                                            label=refs['H1998'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None',
                                            #fillstyle='none',
                                            #mew=2,  # marker edge width
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(CW2004['Age'],
                                            CW2004['Q'],
                                            xerr=CW2004['Age_error'],
                                            yerr=CW2004['Q_error'],
                                            c=colors['yellow'],
                                            label=refs['CW2004'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None',
                                            #fillstyle='none',
                                            #mew=2,  # marker edge width
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(N2008['Age'],
                                            N2008['Q'],
                                            xerr=N2008['Age_error'],
                                            yerr=N2008['Q_error'],
                                            c=colors['navy'],
                                            label=refs['N2008'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None',
                                            #fillstyle='none',
                                            #mew=2,  # marker edge width
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(KH2008['Age'],
                                            KH2008['Q'],
                                            xerr=KH2008['Age_error'],
                                            yerr=KH2008['Q_error'],
                                            c=colors['turquoise'],
                                            label=refs['KH2008'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None',
                                            #fillstyle='none',
                                            #mew=2,  # marker edge width
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(Wright2013['Age'],
                                            Wright2013['Q'],
                                            xerr=Wright2013['Age_error'],
                                            yerr=Wright2013['Q_error'],
                                            c=colors['purple'],
                                            label=refs['Wright2013'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None',
                                            #fillstyle='none',
                                            #mew=2,  # marker edge width
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(Parker2014['Age'],
                                            Parker2014['Q'],
                                            xerr=Parker2014['Age_error'],
                                            yerr=Parker2014['Q_error'],
                                            c=colors['gray'],
                                            label=refs['Parker2014'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None',
                                            #fillstyle='none',
                                            #mew=2,  # marker edge width
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    markers, caps, bars = pyplot.errorbar(Parker2017['Age'],
                                            Parker2017['Q'],
                                            xerr=Parker2017['Age_error'],
                                            yerr=Parker2017['Q_error'],
                                            c=colors['green'],
                                            label=refs['Parker2017'],
                                            marker='D',
                                            ms=12,
                                            elinewidth=2,
                                            capsize=2,
                                            ls='None',
                                            #fillstyle='none',
                                            #mew=2,  # marker edge width
                                            )
    [bar.set_alpha(0.5) for bar in bars]

    for i, txt in enumerate(CW2004['Region']):
        if txt == 'Ophiuchus':
            pyplot.annotate(txt,
                                (CW2004[i]['Age'] + 0.1, CW2004[i]['Q'] - 0.03),
                                fontsize=20,
                                color=colors['yellow'])
        else:
            pyplot.annotate(txt,
                            (CW2004[i]['Age'] + 0.05, CW2004[i]['Q'] + 0.01),
                            fontsize=20,
                            color=colors['yellow'])

    for i, txt in enumerate(Parker2017['Region']):
        if txt == 'NGC1333':
            pyplot.annotate(txt,
                                (Parker2017[i]['Age'] + 0.1, Parker2017[i]['Q'] + 0.01),
                                fontsize=20,
                                color=colors['green'])
        else:
            pyplot.annotate(txt,
                            (Parker2017[i]['Age'] - 0.9, Parker2017[i]['Q'] - 0.03),
                            fontsize=20,
                            color=colors['green'])

    for i, txt in enumerate(Wright2013['Region']):
        pyplot.annotate(txt,
                            (Wright2013[i]['Age'] + 0.1,  Wright2013[i]['Q'] + 0.01),
                            fontsize=20,
                            color=colors['purple'])

    for i, txt in enumerate(Parker2014['Region']):
        pyplot.annotate(txt,
                            (Parker2014[i]['Age'], Parker2014[i]['Q'] + 0.01),
                            fontsize=20,
                            color=colors['gray'])

    for i, txt in enumerate(H1998['Region']):
        pyplot.annotate(txt,
                            (H1998[i]['Age'] - 0.85, H1998[i]['Q'] - 0.01),
                            fontsize=20,
                            color=colors['pink'])

    for i, txt in enumerate(N2008['Region']):
        pyplot.annotate(txt,
                            (N2008[i]['Age'], N2008[i]['Q'] + 0.01),
                            fontsize=20,
                            color=colors['navy'])

    for i, txt in enumerate(KH2008['Region']):
        pyplot.annotate(txt,
                            (KH2008[i]['Age'], KH2008[i]['Q'] + 0.01),
                            fontsize=20,
                            color=colors['turquoise'])

    #pyplot.legend(fontsize=18)
    pyplot.xlabel('Age [Myr]')
    pyplot.ylabel('Q', fontsize=24)
    pyplot.legend(fontsize=18, loc='lower right')

    if save:
        pyplot.savefig('{0}/Q_vs_age.png'.format(save_path))

def main(open_path, nruns, save, save_path):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    #Q_scatter(open_path, nruns, save, save_path)
    Q_vs_time(open_path, nruns, save, save_path)
    #Q_parameter(open_path, nruns)

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

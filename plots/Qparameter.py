from matplotlib import pyplot
import os

from amuse.lab import *
from amuse import io

from mycolors import *
from legends import *


def Q_vs_time(open_path, nruns, save, save_path):
    fig = pyplot.figure()

    times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    Qparam = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

    for n in range(nruns):
        path = '{0}/{1}/disks/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        for f in files[1:]:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            t = float(f.split('t')[1].split('.hdf5')[0])
            q = stars.Qparameter()
            times[n].append(t)
            Qparam[n].append(q)

    for n in range(nruns):
        pyplot.plot(times[n],
                    Qparam[n],
                    lw=3,
                    c=runcolors[n],
                    label=r'Run \#{0}'.format(n))

    pyplot.xlabel('Time [Myr]')
    pyplot.ylabel('Q parameter')

    pyplot.legend()

    if save:
        pyplot.savefig('{0}/Qparam.png'.format(save_path))


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
        converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), 3.0 | units.parsec)
        stars = stars.bound_subset(tidal_radius=Rvir, unit_converter=converter)
        Q = stars.Qparameter()

        print "Run {0}, t = {1:.3f} Myr, Q = {2}".format(n,
                                                          t.value_in(units.Myr),
                                                          Q)


def main(open_path, nruns, save, save_path):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    #Q_vs_time(open_path, nruns, save, save_path)
    Q_parameter(open_path, nruns)


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

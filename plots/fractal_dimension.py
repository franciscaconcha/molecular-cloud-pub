from matplotlib import pyplot
import os

from amuse.lab import *
from amuse import io

from mycolors import *
from legends import *


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


def main(open_path, nruns, save, save_path):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    #fractal_dimension(open_path, nruns, save, save_path)
    fd_vs_time(open_path, nruns, save, save_path)

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

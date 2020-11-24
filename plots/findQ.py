import os

from amuse.lab import *
from amuse import io

from mycolors import *


def Q_vs_time(open_path, nruns, save, save_path):

    times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    Qparam = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

    for n in range(nruns):
        path = '{0}/{1}/disks/'.format(open_path, n)
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        files = [x for x in files if '.hdf5' in x]
        files.sort(key=lambda f: float(filter(str.isdigit, f)))

        print path
        #files = files[0::10]

        for f in files[0::10]:
            stars = read_set_from_file(path + f, 'hdf5', close_file=True)
            t = float(f.split('t')[1].split('.hdf5')[0])
            Rvir = stars.virial_radius()
            converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), Rvir)
            stars = stars.bound_subset(tidal_radius=Rvir, unit_converter=converter)
            q = stars.Qparameter()
            times[n].append(t)
            Qparam[n].append(q)

    print times
    print Qparam

    print "plummer:"

    plummert, plummerq = [], []

    for f in files[1:]:
        stars = read_set_from_file(path + f, 'hdf5', close_file=True)
        t = float(f.split('t')[1].split('.hdf5')[0])
        Rvir = stars.virial_radius()
        converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), Rvir)
        stars = stars.bound_subset(tidal_radius=Rvir, unit_converter=converter)
        q = stars.Qparameter()
        plummert.append(t)
        plummerq.append(q)

    print plummert
    print plummerq



def main(open_path, nruns, save, save_path):

    Q_vs_time(open_path, nruns, save, save_path)




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

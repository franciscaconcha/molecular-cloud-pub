import numpy
import os
from matplotlib import pyplot
from matplotlib.colors import LogNorm

from amuse.lab import *
from amuse.community.fi.interface import *
from prepare_figure import single_frame
from hydrodynamics_class import Hydro

def make_map(sph, N=100, L=1):
    x, y = numpy.indices((N + 1, N + 1))
    x = L * (x.flatten() - N / 2.) / N
    y = L * (y.flatten() - N / 2.) / N
    z = x * 0.
    vx = 0. * x
    vy = 0. * x
    vz = 0. * x

    x = units.parsec(x)
    y = units.parsec(y)
    z = units.parsec(z)
    vx = units.kms(vx)
    vy = units.kms(vy)
    vz = units.kms(vz)

    rho, rhovx, rhovy, rhovz, rhoe = sph.get_hydro_state_at_point(x, y, z, vx, vy, vz)
    rho = rho.reshape((N + 1, N + 1))
    return rho


def plot_molecular_cloud(gasfile, sinksfile, starsfile, path, save_path, prev_file, L=10.0):
    x_label = "x [pc]"
    y_label = "y [pc]"
    pyplot.close('all')
    fig = single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12)

    gas = read_set_from_file(path + gasfile, "amuse")
    sinks = read_set_from_file(path + sinksfile, "hdf5", close_file=True)
    stars = read_set_from_file(path + starsfile, "hdf5", close_file=True)

    print "Ngas={0}, Nsinks={1}, Nstars={2}".format(len(gas), len(sinks), len(stars))

    sph = Hydro(Fi, gas)
    time = str('%.2f' % (sph.gas_particles.get_timestamp().value_in(units.Myr)))

    rho = make_map(sph, N=200, L=L)

    cax = pyplot.imshow(numpy.log10(1.e-5 + numpy.transpose(rho).value_in(units.amu / units.cm ** 3)),
                        extent=[-L / 2, L / 2, -L / 2, L / 2],
                        vmin=1, vmax=5, origin="lower", cmap=pyplot.get_cmap('YlOrBr'))

    """delta = 0.025
    x = y = numpy.arange(-3.0, 3.0, delta)
    X, Y = numpy.meshgrid(x, y)
    Z = numpy.exp(-(X-1) ** 2 - (Y) ** 2)
    #Z2 = numpy.exp(-(X + 2 - 1) ** 2 - (Y - 1) ** 2)
    #Z = (Z1 - Z2) * 2

    cax = pyplot.imshow(Z, interpolation='bilinear', cmap=pyplot.get_cmap('YlOrBr'),
                   origin='upper', extent=[-L / 2, L / 2, -L / 2, L / 2],
                   vmax=abs(Z).max(), vmin=-abs(Z).max())"""


    cbar = fig.colorbar(cax, orientation='vertical', fraction=0.045)
    cbar.set_label('projected density [$amu/cm^3$]', rotation=270, labelpad=25)

    if len(sinks):
        pyplot.scatter(sinks.x.value_in(units.parsec),
                       sinks.y.value_in(units.parsec),
                       c='white', alpha=0.5)

    if len(stars):
        pyplot.scatter(stars.x.value_in(units.parsec),
                       stars.y.value_in(units.parsec),
                       c='blue', marker="*", alpha=0.5)

    pyplot.xlim(-L / 2., L / 2.)
    pyplot.ylim(-L / 2., L / 2.)
    pyplot.title("Molecular cloud at time={0} Myr".format(time))
    pyplot.xlabel("x [pc]")
    pyplot.ylabel("y [pc]")
    ff = gasfile.split('_')[-1]
    pyplot.savefig('{0}/{1}.png'.format(save_path, ff[1:].split('.')[0]))
    #pyplot.show()


def main(filename, path, save_path):
    prev_file = None

    if path is not None:
        files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'

        gas_files = [x for x in files if 'gas' in x]
        gas_files.sort(key=lambda f: int(filter(str.isdigit, f)))

        sink_files = [x for x in files if 'gravity_sinks' in x]
        sink_files.sort(key=lambda f: int(filter(str.isdigit, f)))

        star_files = [x for x in files if 'gravity_stars' in x]
        star_files.sort(key=lambda f: int(filter(str.isdigit, f)))

        for g, si, ss in zip(gas_files, sink_files, star_files):
            filename = '{0}/{1}'.format(path, g)
            plot_molecular_cloud(g, si, ss, path, save_path, prev_file)
            prev_file = filename

    else:
        plot_molecular_cloud(filename, save_path, prev_file)


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default="GMC_R2pcN20k_SE_T45Myr.amuse",
                      help="input filename [%default]")
    result.add_option("-p", dest="path", default=None,
                      help="path for input files [%default]")
    result.add_option("-s", dest="save_path", default=".",
                      help="save path [%default]")
    return result


if __name__ in ('__main__', '__plot__'):
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)






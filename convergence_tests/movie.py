import numpy
import os.path
from matplotlib import pyplot

from amuse.lab import *
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


def plot_molecular_cloud(filename, save_path, L=10.0):
    x_label = "x [pc]"
    y_label = "y [pc]"
    fig = single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12)

    print "read file:", filename
    gas = read_set_from_file(filename, "amuse")
    sinkfile = filename.split("_gas_")[0] + "_sink_" + filename.split("_gas_")[1]
    starfile = filename.split("_gas_")[0] + "_stars_" + filename.split("_gas_")[1]
    print sinkfile
    print starfile

    if os.path.isfile(sinkfile):
        sinks = read_set_from_file(sinkfile, "hdf5", close_file=True)
    else:
        sinks = Particles(0)

    if os.path.isfile(starfile):
        stars = read_set_from_file(starfile, "hdf5", close_file=True)
    else:
        stars = Particles(0)

    print "Ngas={0}, Nsinks={1}, Nstars={2}".format(len(gas), len(sinks), len(stars))

    sph = Hydro(Fi, gas)
    time = 0 | units.Myr

    rho = make_map(sph, N=200, L=L)

    cax = pyplot.imshow(numpy.log10(1.e-5 + rho.value_in(units.amu / units.cm ** 3)),
                        extent=[-L / 2, L / 2, -L / 2, L / 2],
                        vmin=1, vmax=5, origin="lower", cmap=pyplot.get_cmap('YlOrBr'))

    cbar = fig.colorbar(cax, orientation='vertical', fraction=0.045)
    cbar.set_label('projected density [$amu/cm^3$]', rotation=270, labelpad=25)

    cm = pyplot.cm.get_cmap('Greys')
    if len(sinks):
        m = 100 * sinks.radius.value_in(units.parsec)
        #100 * numpy.log10(sinks.mass / sinks.mass.min())
        c = numpy.sqrt(sinks.mass / sinks.mass.max())
        pyplot.scatter(sinks.y.value_in(units.parsec), sinks.x.value_in(units.parsec),
                       c=c, s=m, lw=0, cmap=cm)

    cm = pyplot.cm.get_cmap('cool')
    if len(stars):
        m = 100 * numpy.log10(stars.mass / stars.mass.min())
        c = numpy.sqrt(stars.mass / stars.mass.max())
        pyplot.scatter(stars.y.value_in(units.parsec), stars.x.value_in(units.parsec),
                       c=c, marker="*", alpha=0.5, lw=0, cmap=cm)

    pyplot.xlim(-L / 2., L / 2.)
    pyplot.ylim(-L / 2., L / 2.)
    pyplot.title("Molecular cloud at time=" + time.as_string_in(units.Myr))
    pyplot.xlabel("x [pc]")
    pyplot.ylabel("x [pc]")
    pyplot.title("GMC at time=" + time.as_string_in(units.Myr))
    ff = filename.split('_')[-1]
    #pyplot.savefig('{0}/{1}.png'.format(save_path, ff[1:].split('.')[0]))
    pyplot.show()


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default="GMC_R2pcN20k_SE_T45Myr.amuse",
                      help="input filename [%default]")
    result.add_option("-s", dest="save_path", default=".",
                      help="save path [%default]")
    return result


if __name__ in ('__main__', '__plot__'):
    o, arguments = new_option_parser().parse_args()
    plot_molecular_cloud(**o.__dict__)






import numpy
import os.path
from matplotlib import pyplot

from amuse.lab import *
from prepare_figure import single_frame
from hydrodynamics_class import Hydro


def make_map(sph, N=100, L=1):
    x, y = numpy.indices((N+1, N+1))
    x = L * (x.flatten() - N/2.) / N
    y = L * (y.flatten() - N/2.) / N
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
    rho = rho.reshape((N+1, N+1))
    return rho


def plot_molecular_cloud(filename, save_path):
    print "read file:", filename
    gas = read_set_from_file(filename, "amuse")
    sinkfile = filename.split("_gas_")[0] + "_sink_" + filename.split("_gas_")[1]
    starfile = filename.split("_star_")[0] + "_sink_" + filename.split("_gas_")[1]

    if os.path.isfile(sinkfile):
        sinks = read_set_from_file(sinkfile, "hdf5", close_file=True)
    else:
        sinks = Particles(0)

    if os.path.isfile(starfile):
        stars = read_set_from_file(starfile, "hdf5", close_file=True)
    else:
        stars = Particles(0)

    print "N=", len(gas), len(sinks), len(stars)
    #print sinks.mass.in_(units.MSun)
    
    #sinks =  bodies.history.next()
    L = 10.0
#    if len(sinks):
#        sinks = sinks.select(lambda r: r.length()<0.5*L|units.parsec,["position"])
    #print sinks
    hydro = Hydro(Fi, gas)
    time = 0 | units.Myr

    if len(sinks) > 0:
        hydro.code.dm_particles.add_particles(sinks)

    if len(stars) > 0:
        hydro.code.dm_particles.add_particles(sinks)

    #plot_hydro_and_stars(filename, save_path, time, hydro, L=10)
    plot_hydro(filename, save_path, time, hydro, L)


def make_stars(cluster_particle):
    sfe = 0.3
    mmean = 1.0 | units.MSun
    N = int(sfe*cluster_particle.mass/mmean)
    stars = Particles(0)
    print "N_cluster=", N
    if N>0:
        masses = new_salpeter_mass_distribution(N, 0.3|units.MSun, min(100|units.MSun, cluster_particle.mass))

        r = cluster_particle.h_smooth  
        converter=nbody_system.nbody_to_si(masses.sum(),r)
        stars = new_plummer_model(N, convert_nbody=converter)
        stars.mass = masses
        stars.position += cluster_particle.position
        stars.velocity += cluster_particle.velocity
    return stars


def get_stars_from_molecular_cloud(parts):
    cutoff_density = 10000 | units.amu/units.cm**3
    stars = Particles(0)
    for ip in parts:
        if ip.rho>cutoff_density:
            local_stars = make_stars(ip)
            stars.add_particles(local_stars)
    return stars


def plot_hydro_and_stars(filename, save_path, time, sph, L=10):
    fig = pyplot.figure(figsize=(12, 12))
    rho = make_map(sph, N=200, L=L)
    cax = pyplot.imshow(numpy.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)),
                        extent=[-L/2, L/2, -L/2, L/2], vmin=1, vmax=5, origin="lower")
    cbar = fig.colorbar(cax, ticks=[4, 7.5, 11], orientation='vertical', fraction=0.045)
    cbar.set_label('projected density [$amu/cm^3$]', rotation=270, labelpad=25)
    pyplot.tight_layout()

    #stars = get_stars_from_molecular_cloud(sph.gas_particles)
    stars = sph.dm_particles

    if len(stars):
        #m =  100.0*stars.mass/max(stars.mass)
        m = 100.0 * stars.mass / stars.mass.max()
        c = 'k' #stars.mass/stars.mass.mean()
        x = -stars.x.value_in(units.parsec)
        y = stars.y.value_in(units.parsec)
        pyplot.scatter(x, y, s=m, c=c, lw=0)

    pyplot.xlim(-L/2., L/2.)
    pyplot.ylim(-L/2., L/2.)
    pyplot.title("Molecular cloud at time="+time.as_string_in(units.Myr))
    pyplot.xlabel("x [pc]")
    pyplot.ylabel("x [pc]")
    pyplot.title("GMC at time="+time.as_string_in(units.Myr))

    #pyplot.show()
    ff = filename.split('_')[-1]
    pyplot.savefig('{0}/{1}.png'.format(save_path, ff[1:].split('.')[0]))


def plot_hydro(filename, save_path, time, sph, L=10):
    x_label = "x [pc]"
    y_label = "y [pc]"
    fig = single_frame(x_label, y_label, logx=False, logy=False, xsize=12, ysize=12)

    gas = sph.code.gas_particles
    dmp = sph.code.dm_particles
    rho=make_map(sph, N=200, L=L)

    cax = pyplot.imshow(numpy.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)),
                        extent=[-L/2, L/2, -L/2, L/2], vmin=1, vmax=5, origin="lower")

    cbar = fig.colorbar(cax, orientation='vertical', fraction=0.045)
    cbar.set_label('projected density [$amu/cm^3$]', rotation=270, labelpad=25)

    cm = pyplot.cm.get_cmap('RdBu')
#    cm = pyplot.cm.jet #gist_ncar

    if len(dmp):
        m = 100 * numpy.log10(dmp.mass / dmp.mass.min())
        c = numpy.sqrt(dmp.mass / dmp.mass.max())
        pyplot.scatter(dmp.y.value_in(units.parsec), dmp.x.value_in(units.parsec), c=c, s=m, lw=0, cmap=cm)

    #pyplot.show()
    ff = filename.split('_')[-1]
    pyplot.savefig('{0}/{1}.png'.format(save_path, ff[1:].split('.')[0]))


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



    


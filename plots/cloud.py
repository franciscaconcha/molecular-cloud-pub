import numpy
from matplotlib import pyplot
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid



from amuse.community.fi.interface import FiMap
from amuse.io import read_set_from_file
from amuse.units import units, nbody_system, constants


mH = 1.008 * constants.u


def column_density (gas_particles, stars, ax, N=480, box_size=10.|units.parsec):
    '''
    Make a column density plot of a distribution of AMUSE gas particles.
    Column density is volume density integrated along a line.

    gas_particles: gas particle set to make plot of, must at least have mass attribute
        (amuse particle set)
    N: number of pixels on each side (int)
    box_size: physical size of plotted region (scalar, unit length)
    '''

    converter = nbody_system.nbody_to_si(gas_particles.mass.sum(),
        gas_particles.position.lengths().max())

    mapper = FiMap(converter)

    # Compute column density over parallell lines-of-sight; effectively look from infinitely far
    mapper.parameters.projection_mode = 'parallel'

    # Physical size of imaged region
    mapper.parameters.image_width = box_size

    # Number of image pixels
    mapper.parameters.image_size = (N, N)

    # Coordinates of image center
    mapper.parameters.image_target = [0., 0., 0.] | units.parsec

    # Camera position (unit vector)
    mapper.parameters.projection_direction = [0, 0, -1]

    # Image orientation
    mapper.parameters.upvector = [0, 1, 0]

    # Particles to include in image
    mapper.particles.add_particles(gas_particles)

    # Quantity to make image of
    weight = gas_particles.mass.value_in(units.MSun)

    # Make image; every pixel contains the total mass within that pixel
    image = mapper.image.pixel_value.T

    mapper.stop()

    # Dividing by pixel size gives mass column density
    bs = box_size.value_in(units.parsec)
    surface_density = image / (bs/N)**2

    if numpy.max(surface_density) == 0.:
        print ("[WARNING] image identically 0")


    p = ax.imshow(surface_density, norm=LogNorm(),
        extent=[-bs/2., bs/2., -bs/2., bs/2.], origin='lower', cmap='RdPu',
        #vmin=numpy.max(surface_density)/1e3, vmax=numpy.max(surface_density))
        vmin=100, vmax=1E4)
    # Since the color scale is logarithmic and calibrated on the max, this does not work if
    # the image is identically zero!

    if len(stars):
        ax.scatter(stars.x.value_in(units.parsec),
                       stars.y.value_in(units.parsec),
                       s=10,
                       marker=".",
                       alpha=0.2,
                       c='black',
                       #facecolors='none',
                       )

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')

    ax.set_title('t = {0:.2f} Myr'.format(gas_particles.get_timestamp().value_in(units.Myr)))
    return p


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-p", dest="path", default=".",
                      help="path for files [%default]")
    result.add_option("-S", dest="save", type="int", default=0,
                      help="save plot? [%default]")
    result.add_option("-s", dest="save_path", type="string", default='.',
                      help="path to save the results [%default]")
    return result


if __name__ in ('__main__', '__plot__'):
    o, arguments = new_option_parser().parse_args()

    pyplot.style.use('paper')

    indices = [1, 16, 18]

    fig, axs = pyplot.subplots(nrows=1,
                               ncols=4,
                               #sharey=True,
                               figsize=(18, 8),
                               gridspec_kw={"width_ratios": [1, 1, 1, 0.1]})
    fig.subplots_adjust(wspace=0.4)

    j = 0

    for fileindex in indices:
        print fileindex
        ax = axs[j]
        gas_particles = read_set_from_file(
            o.path + 'hydro_gas_particles_i{:04}.amuse'.format(int(fileindex)),
            'hdf5',
            close_file=True)
        print gas_particles.get_timestamp().value_in(units.Myr)

        try:
            stars = read_set_from_file(
                o.path +'hydro_stars_particles_i{:04}.amuse'.format(int(fileindex)), 'hdf5')
            print stars.get_timestamp().value_in(units.Myr)

        except:
            stars = []

        im = column_density(gas_particles, stars, ax)
        j += 1

    cbar = fig.colorbar(im, cax=axs[j], fraction=0.5, pad=0.04)
    cbar.set_label('Column Density [M$_\\odot$ pc$^{-2}$]')

    pyplot.tight_layout()

    if o.save:
        pyplot.savefig('{0}/cloud.png'.format(o.save_path))
    else:
        pyplot.show()

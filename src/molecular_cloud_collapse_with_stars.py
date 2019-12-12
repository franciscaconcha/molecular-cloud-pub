import numpy
from scipy import interpolate

from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from amuse.couple.bridge import Bridge
from amuse.ic.fractalcluster import new_fractal_cluster_model

from cooling_class import SimplifiedThermalModel, SimplifiedThermalModelEvolver
from hydrodynamics_class import Hydro
from gravity_class import Gravity

######## FRIED grid ########
# Yes doing this with global variables is bad practice... but practical since
# these values will be immutable and I will use them 100s of times

# Read FRIED grid
grid = numpy.loadtxt('data/friedgrid.dat', skiprows=2)

# Getting only the useful parameters from the grid (not including Mdot)
FRIED_grid = grid[:, [0, 1, 2, 4]]
grid_log10Mdot = grid[:, 5]

grid_stellar_mass = FRIED_grid[:, 0]
grid_FUV = FRIED_grid[:, 1]
grid_disk_mass = FRIED_grid[:, 2]
grid_disk_radius = FRIED_grid[:, 3]


def write_data(path, timestamp, hydro, index=0, stars=Particles(0)):
    hydro.write_set_to_file(path, index=index)
    filename = "{0}/hydro_stars_particles_i{1:04}.amuse".format(path, index)
    if len(stars) > 0:
        write_set_to_file(stars, filename, "hdf5",
                          timestamp=timestamp,
                          append_to_file=False)


def generate_initial_conditions_for_molecular_cloud(N, Mcloud, Rcloud):
    conv = nbody_system.nbody_to_si(Mcloud, Rcloud)
    gas = molecular_cloud(targetN=N, convert_nbody=conv,
                          base_grid=body_centered_grid_unit_cube).result
    gas.name = "gas"

    return gas


def make_stars_from_sink(sink, stellar_mass, time):

    # Delay time for next star formation is a decay from the tff of the sink
    delay_time = sink.tff * numpy.exp(-0.1 * time.value_in(units.Myr))

    print "Forming star of mass {0} from sink mass {1}".format(stellar_mass.in_(units.MSun),
                                                               sink.mass.in_(units.MSun))
    # If sink is massive enough and it's time to form a star
    stars_from_sink = Particles(1)
    stars_from_sink.stellar_mass = stellar_mass
    stars_from_sink.disk_mass = 0.1 * stars_from_sink.stellar_mass
    stars_from_sink.mass = stars_from_sink.stellar_mass + stars_from_sink.disk_mass
    sink.mass -= stars_from_sink.mass

    # TODO create disks here, or ourside this function? disk_codes could be global

    # Find position offset inside sink radius
    Rsink = sink.radius.value_in(units.parsec)
    offset = numpy.random.uniform(-Rsink, Rsink) | units.parsec
    stars_from_sink.x = sink.x + offset
    stars_from_sink.y = sink.y + offset
    stars_from_sink.z = sink.z + offset

    stars_from_sink.vx = sink.vx
    stars_from_sink.vy = sink.vy
    stars_from_sink.vz = sink.vz

    return stars_from_sink, delay_time


def find_indices(column,
                 val):
    """
    Return indices of column values in between which val is located.
    Return i,j such that column[i] < val < column[j]

    :param column: column where val is to be located
    :param val: number to be located in column
    :return: i, j indices
    """

    # The largest element of column less than val
    try:
        value_below = column[column < val].max()
    except ValueError:
        # If there are no values less than val in column, return smallest element of column
        value_below = column.min()
    # Find index
    index_i = numpy.where(column == value_below)[0][0]

    # The smallest element of column greater than val
    try:
        value_above = column[column > val].min()
    except ValueError:
        # If there are no values larger than val in column, return largest element of column
        value_above = column.max()
    # Find index
    index_j = numpy.where(column == value_above)[0][0]

    return int(index_i), int(index_j)


def get_disk_radius(disk,
                    density_limit=1E-10):
    """ Calculate the radius of a disk in a vader grid.

    :param disk: vader disk
    :param density_limit: density limit to designate disk border
    :return: disk radius in units.au
    """
    prev_r = disk.grid[0].r

    for i in range(len(disk.grid.r)):
        cell_density = disk.grid[i].column_density.value_in(units.g / units.cm ** 2)
        if cell_density < density_limit:
            return prev_r.value_in(units.au) | units.au
        prev_r = disk.grid[i].r

    return prev_r.value_in(units.au) | units.au


def get_disk_mass(disk,
                  radius):
    """ Calculate the mass of a vader disk inside a certain radius.

    :param disk: vader disk
    :param radius: disk radius to consider for mass calculation
    :return: disk mass in units.MJupiter
    """
    mass_cells = disk.grid.r[disk.grid.r <= radius]
    total_mass = 0

    for m, d, a in zip(mass_cells, disk.grid.column_density, disk.grid.area):
        total_mass += d.value_in(units.MJupiter / units.cm**2) * a.value_in(units.cm**2)

    return total_mass | units.MJupiter


def get_disk_density(disk):
    """ Calculate the mean density of the disk, not considering the outer, low density limit.

    :param disk: vader disk
    :return: mean disk density in g / cm**2
    """
    radius = get_disk_radius(disk)
    radius_index = numpy.where(disk.grid.r.value_in(units.au) == radius.value_in(units.au))
    density = disk.grid[:radius_index[0][0]].column_density.value_in(units.g / units.cm**2)
    return numpy.mean(density) | (units.g / units.cm**2)


# TODO write separate function to calculate total radiation over star
def photoevaporation_mass_loss(star, radiation, dt):
    global FRIED_grid, grid_log10Mdot, grid_stellar_mass, grid_FUV, grid_disk_mass, grid_disk_radius

    if star.EUV:
        # Photoevaporative mass loss in MSun/yr. Eq 20 from Johnstone, Hollenbach, & Bally 1998
        # From the paper: e ~ 3, x ~ 1.5
        photoevap_Mdot = 2. * 1E-9 * 3 * 4.12 * (star.disk_radius.value_in(units.cm) / 1E14)

        # Calculate total mass lost due to EUV photoevaporation during dt, in MSun
        total_photoevap_mass_loss_euv = float(photoevap_Mdot * dt.value_in(units.yr)) | units.MSun

        # Back to False for next time
        star.EUV = False  # TODO after calling photoevaporation_mass_loss, stars' EUV should be set to False
    else:
        total_photoevap_mass_loss_euv = 0.0 | units.MSun

    # FUV regime -- Use FRIED grid

    # For the small star, I want to interpolate the photoevaporation mass loss
    # xi will be the point used for the interpolation. Adding star values...
    xi = numpy.ndarray(shape=(1, 4), dtype=float)
    xi[0][0] = star.stellar_mass.value_in(units.MSun)
    xi[0][1] = radiation
    # TODO should disk_codes also be global in this script?
    xi[0][3] = get_disk_radius(disk_codes[disk_codes_indices[star.key]]).value_in(units.au)
    xi[0][2] = get_disk_mass(disk_codes[disk_codes_indices[star.key]], xi[0][3] | units.au).value_in(units.MJupiter)

    # Building the subgrid (of FRIED grid) over which I will perform the interpolation
    subgrid = numpy.ndarray(shape=(8, 4), dtype=float)

    # Finding indices between which star.mass is located in the grid
    stellar_mass_i, stellar_mass_j = find_indices(grid_stellar_mass, star.stellar_mass.value_in(units.MSun))

    subgrid[0] = FRIED_grid[stellar_mass_i]
    subgrid[1] = FRIED_grid[stellar_mass_j]

    # Finding indices between which the radiation over the small star is located in the grid
    FUV_i, FUV_j = find_indices(grid_FUV, total_radiation[star.key])
    subgrid[2] = FRIED_grid[FUV_i]
    subgrid[3] = FRIED_grid[FUV_j]

    # Finding indices between which star.disk_mass is located in the grid
    disk_mass_i, disk_mass_j = find_indices(grid_disk_mass, star.disk_mass.value_in(units.MJupiter))
    subgrid[4] = FRIED_grid[disk_mass_i]
    subgrid[5] = FRIED_grid[disk_mass_j]

    # Finding indices between which star.disk_radius is located in the grid
    disk_radius_i, disk_radius_j = find_indices(grid_disk_radius, star.disk_radius.value_in(units.au))
    subgrid[6] = FRIED_grid[disk_radius_i]
    subgrid[7] = FRIED_grid[disk_radius_j]

    # Adding known values of Mdot, in the indices found above, to perform interpolation
    Mdot_values = numpy.ndarray(shape=(8,), dtype=float)
    indices_list = [stellar_mass_i, stellar_mass_j,
                    FUV_i, FUV_j,
                    disk_mass_i, disk_mass_j,
                    disk_radius_i, disk_radius_j]
    for x in indices_list:
        Mdot_values[indices_list.index(x)] = grid_log10Mdot[x]

    # Interpolate!
    # Photoevaporative mass loss in log10(MSun/yr)
    photoevap_Mdot = interpolate.griddata(subgrid, Mdot_values, xi, method="nearest")

    # Calculate total mass lost due to photoevaporation during dt, in MSun
    total_photoevap_mass_loss_fuv = float(numpy.power(10, photoevap_Mdot) * dt.value_in(units.yr)) | units.MSun

    return total_photoevap_mass_loss_euv + total_photoevap_mass_loss_fuv


def run_molecular_cloud(gas_particles, sink_particles, SFE, method, tstart, tend, dt_diag, save_path, index=0):

    hydro = Hydro(Fi, gas_particles, tstart)

    if len(sink_particles) > 0:
        hydro.sink_particles.add_particles(sink_particles)
        hydro.sink_particles.synchronize_to(hydro.code.dm_particles)

    gravity = None
    gravhydro = None
    gravity_sinks = None

    dt = min(dt_diag, 0.1 | units.Myr)
    t_diag = 0 | units.Myr

    E0 = hydro.gas_particles.kinetic_energy() \
         + hydro.gas_particles.potential_energy() \
         + hydro.gas_particles.thermal_energy()
    time = gas_particles.get_timestamp()

    # Sample IMF for single star formation
    IMF_masses = new_kroupa_mass_distribution(10000, mass_max=50 | units.MSun)
    current_mass = 0  # To keep track of formed stars

    stars = Particles(0)

    sink_formation = True  # To keep track of SFE
    local_sinks = Particles(0)

    # Channels from hydro to local sinks and viceversa
    hydro_sinks_to_framework = hydro.sink_particles.new_channel_to(local_sinks)
    framework_to_hydro_sinks = local_sinks.new_channel_to(hydro.sink_particles)

    # To keep track of low mass (disked) and high mass (radiating) stars
    low_mass, high_mass = [], []

    while time < tend:
        time += dt
        print "Evolve to time=", time.in_(units.Myr)
        Mtot = 0 | units.MSun

        if sink_formation:  # This means that SPH code is still active. Might rename this variable later.
            if gravity is None:
                if len(hydro.sink_particles) > 0:
                    Mtot = hydro.gas_particles.mass.sum() + hydro.sink_particles.mass.sum()
                else:
                    Mtot = hydro.gas_particles.mass.sum()
            else:  # If there's a gravity code running, we count the mass of the stars as well
                # TODO fix this for when gravity_sinks is running
                Mtot = hydro.gas_particles.mass.sum() + hydro.sink_particles.mass.sum() + gravity.particles.mass.sum()

            if len(local_sinks) > 0:
                MC_SFE = local_sinks.mass.sum() / Mtot
            else:
                MC_SFE = 0

            if MC_SFE >= SFE:
                print "************** SFE reached, sinks will stop forming **************"
                sink_formation = False
                # TODO stop hydro code, kick out all gas, keep going with Nbody

                # Final synchro and copy of hydro sinks to local sinks
                hydro.sink_particles.synchronize_to(local_sinks)
                hydro_sinks_to_framework.copy()

                # Create new gravity code for sinks only # TODO check: is this a good idea?
                gravity_sinks = Gravity(ph4, local_sinks)
                gravity_sinks_to_framework = gravity_sinks.code.particles.new_channel_to(local_sinks)
                framework_to_gravity_sinks = local_sinks.new_channel_to(gravity_sinks.code.particles)

                # From this point on I will simply keep evolving the gravity and gravity_sinks codes

        for sink in local_sinks:
            # Iterate over local_sinks instead of hydro.sink_particles so that the leftover
            # sinks can still form stars after the gas code is stopped.

            if sink.mass > IMF_masses[current_mass] and sink.form_star:
                stars_from_sink, delay_time = make_stars_from_sink(sink, IMF_masses[current_mass], time)
                sink.form_star = False
                sink.time_threshold = time + delay_time  # Next time at which this sink should form a star

                # Update sink masses and delay times from framework to codes
                if gravity_sinks is None:
                    framework_to_hydro_sinks.copy()
                else:
                    framework_to_gravity_sinks.copy()

                current_mass += 1
                stars.add_particles(stars_from_sink)

                # TODO add disk parameters to stars_from_sink (do it in make_stars_from_sink)
                # TODO create disk codes... here or in make_stars_from_sink?

                # To keep track of "disked" and "FUV-radiating" stars
                if stars_from_sink.stellar_mass > 1.9 | units.MSun:
                    high_mass.append(stars_from_sink.key)
                else:
                    low_mass.append(stars_from_sink.key)  # TODO ojo: Mstar < 0.05 MSun will not have disks


                # I don't care about gravity_sinks for this block
                if gravity is None:
                    print 'Starting gravity code and Bridge'
                    gravity_offset_time = time
                    gravity = Gravity(ph4, stars)
                    gravity_to_framework = gravity.code.particles.new_channel_to(stars)
                    framework_to_gravity = stars.new_channel_to(gravity.code.particles)
                    gravity_to_framework.copy()

                    gravhydro = Bridge()
                    gravhydro.add_system(gravity, (hydro,))
                    gravhydro.add_system(hydro, (gravity,))
                    gravhydro.timestep = 0.1 * dt
                else:
                    gravity.code.particles.add_particles(stars_from_sink)
                    gravity_to_framework.copy()

            elif sink.mass > IMF_masses[current_mass] and not sink.form_star:
                print "Sink is massive enough, but it's not yet time to form a star."
                if time >= sink.time_threshold:
                    sink.form_star = True

                    # Update sink.form_star in corresponding code
                    if gravity_sinks is None:
                        framework_to_hydro_sinks.copy()
                    else:
                        framework_to_gravity_sinks.copy()
                    # sink will form a star in the next timestep

            elif sink.mass < IMF_masses[current_mass] and sink.form_star:
                print "Sink is not massive enough to form this star."
                # sink.form_star = False

            gravity_to_framework.copy()
            if gravity_sinks is None:
                framework_to_hydro_sinks.copy()
            else:
                framework_to_gravity_sinks.copy()

        if gravhydro is None:
            print "Evolving hydro only"
            hydro.evolve_model(time)
            # Synchronize sinks, then update local sinks
            hydro.sink_particles.synchronize_to(local_sinks)
            hydro_sinks_to_framework.copy()
        else:
            if sink_formation:
                print "EVOLVING GRAVHYDRO with {0} particles".format(len(gravity.particles))
                gravhydro.evolve_model(time - gravity_offset_time)
                # Synchronize sinks, then update local sinks
                hydro.sink_particles.synchronize_to(local_sinks)
                gravity.particles.synchronize_to(stars)
                hydro_sinks_to_framework.copy()
                gravity_to_framework.copy()
                print "GRAVHYDRO.MODEL_TIME: {0}".format(gravhydro.model_time.in_(units.Myr))
            else:
                print "EVOLVING GRAVITY AND GRAVITY_SINKS ONLY"
                gravity.evolve_model(time - gravity_offset_time)
                gravity_sinks.evolve_model(time)
                gravity_to_framework.copy()
                gravity_sinks_to_framework.copy()
                #local_sinks.synchronize_to(hydro.sink_particles)

            # TODO also here add stellar evolution, disk evolution, etc

        E = hydro.gas_particles.kinetic_energy() \
            + hydro.gas_particles.potential_energy() \
            + hydro.gas_particles.thermal_energy()
        E_th = hydro.gas_particles.thermal_energy()
        Eerr = (E - E0) / E0
        print 'energy=', E, 'energy_error=', Eerr, 'e_th=', E_th
        print "maximal_density:", gas_particles.rho.max().in_(units.MSun / units.parsec ** 3)

        #hydro.print_diagnostics()
        if gravhydro is None:
            print "No gravhydro yet."
        else:
            print "gravhydro"
            # print_diagnostics(gravhydro)
        if time > t_diag:
            index += 1
            t_diag += dt

            if sink_formation:
                write_data(save_path, time, hydro=hydro, index=index, stars=stars)

            # Saving stars and local sinks separately from the Hydro files
            # Saving star particles
            write_set_to_file(stars,
                              '{0}/gravity_stars_t{1:.2f}Myr.hdf5'.format(save_path,
                                                                          time.value_in(units.Myr)),
                              'hdf5')
            # Saving local sink particles
            write_set_to_file(local_sinks,
                              '{0}/gravity_sinks_t{1:.2f}Myr.hdf5'.format(save_path,
                                                                          time.value_in(units.Myr)),
                              'hdf5')

    #print len(gravity.code.particles)
    hydro.stop(), gravity.stop(), gravity_sinks.stop()
    return gas_particles


def main(filename, save_path, tend, dt_diag, Ncloud, Mcloud, Rcloud, method):
    if len(filename) == 0:
        try:
            import os
            os.makedirs(save_path)
        except OSError, e:
            if e.errno != 17:
                raise
            # time.sleep might help here
            pass

        gas_particles = generate_initial_conditions_for_molecular_cloud(o.Ncloud,
                                                                        Mcloud=o.Mcloud,
                                                                        Rcloud=o.Rcloud)
        rho_cloud = o.Mcloud / o.Rcloud ** 3
        tff = 0.5427 / numpy.sqrt(constants.G * rho_cloud)
        print "Freefall timescale=", tff.in_(units.Myr)
        rho_cloud = 3. * o.Mcloud / (4. * numpy.pi * o.Rcloud ** 3)
        print rho_cloud

        #tend = 10 * tff
        dt_diag = 0.1 * tff
        hydro = Hydro(Fi, gas_particles)
        write_data(save_path, hydro.model_time, hydro)
        filename = "hydro_gas_particles_i{0:04}.amuse".format(0)
        hydro.stop()

    gas_particles = read_set_from_file("{0}/{1}".format(save_path, filename), "hdf5", close_file=True)
    start_time = gas_particles.get_timestamp()

    index = int(filename.split("_i")[1].split(".amuse")[0])
    sinkfile = filename.split("_gas_")[0] + "_sink_" + filename.split("_gas_")[1]

    import os.path
    if os.path.isfile(sinkfile):
        sink_particles = read_set_from_file(sinkfile, "hdf5", close_file=True)
    else:
        sink_particles = Particles(0)

    print "Time= {0}".format(start_time.in_(units.Myr))
    print "index = {0}, Ngas = {1}, Nsinks = {2}".format(index, len(gas_particles), len(sink_particles))

    SFE = 0.4
    #method = 'cluster'
    #method = 'single'

    parts = run_molecular_cloud(gas_particles, sink_particles, SFE, method, start_time, tend, dt_diag, save_path, index)


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename",
                      default="",
                      help="input filename")
    result.add_option("-s", dest="save_path",
                      default="./results/",
                      help="save path for results")
    result.add_option("--tend", dest="tend",
                      unit=units.Myr,
                      type="float",
                      default=5.0 | units.Myr,
                      help="end time")
    result.add_option("--dt_diag", dest="dt_diag",
                      unit=units.Myr,
                      type="float",
                      default=0.1 | units.Myr,
                      help="diagnosticstime step")
    result.add_option("--Ncloud", dest="Ncloud",
                      default=1000,
                      type="float",
                      help="number of gas particles.")
    result.add_option("--Mcloud", dest="Mcloud",
                      unit=units.MSun,
                      type="float",
                      default=1000 | units.MSun,
                      help="cloud mass")
    result.add_option("--Rcloud", dest="Rcloud",
                      unit=units.parsec,
                      type="float",
                      default=3 | units.parsec,
                      help="cloud size")
    result.add_option("-m", dest="method",
                      type="string",
                      default='cluster',
                      help="method for star formation, single or cluster")

    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    numpy.random.seed(3141)
    main(**o.__dict__)

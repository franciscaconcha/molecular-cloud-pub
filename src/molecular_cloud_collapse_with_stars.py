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
from disks_class import Disk


######## FRIED grid ########
# Yes doing this with global variables is bad practice... but practical since
# these values will be immutable and I will use them 100s of times

# Read FRIED grid
grid_data = numpy.loadtxt('data/friedgrid.dat', skiprows=2)

# Getting only the useful parameters from the grid (not including Mdot)
FRIED_grid = grid_data[:, [0, 1, 2, 4]]
grid_log10Mdot = grid_data[:, 5]

grid_stellar_mass = FRIED_grid[:, 0]
grid_FUV = FRIED_grid[:, 1]
grid_disk_mass = FRIED_grid[:, 2]
grid_disk_radius = FRIED_grid[:, 3]

diverged_disks = {}
disk_codes_indices = {}
disk_codes = []


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


def accretion_rate(mass):

    return numpy.power(10, (1.89 * numpy.log10(mass.value_in(units.MSun)) - 8.35)) | units.MSun / units.yr


def make_star_from_sink(sink, stellar_mass, time, alpha=1E-4, ncells=50):

    # Delay time for next star formation is a decay from the tff of the sink
    delay_time = sink.tff * numpy.exp(-0.1 * time.value_in(units.Myr))

    print "Forming star of mass {0} from sink mass {1}".format(stellar_mass.in_(units.MSun),
                                                               sink.mass.in_(units.MSun))
    # If sink is massive enough and it's time to form a star
    new_star = Particles(1)
    new_star.stellar_mass = stellar_mass
    new_star.disk_mass = 0.1 * new_star.stellar_mass
    new_star.mass = new_star.stellar_mass + new_star.disk_mass
    sink.mass -= new_star.mass

    # 'Normal' star parameters: location, velocity, etc
    # Find position offset inside sink radius
    Rsink = sink.radius.value_in(units.parsec)
    offset = numpy.random.uniform(-Rsink, Rsink) | units.parsec
    new_star.x = sink.x + offset
    new_star.y = sink.y + offset
    new_star.z = sink.z + offset

    new_star.vx = sink.vx
    new_star.vy = sink.vy
    new_star.vz = sink.vz

    # Star parameters needed for photoevaporation routines
    new_star.bright = new_star.stellar_mass > 1.9 | units.MSun

    if not new_star.bright:
        new_star.disk_radius = 100 * (new_star.stellar_mass.value_in(units.MSun) ** 0.5) | units.au
        new_star.g0 = 0.0

        new_star.code = True

        # Creating disk for new star
        global disk_codes, disk_codes_indices, diverged_disks
        new_disk = Disk(new_star.disk_radius, new_star.disk_mass, alpha, n_cells=ncells, linear=False)
        s_code = new_disk.code
        print s_code

        s_code.parameters.inner_pressure_boundary_mass_flux = accretion_rate(new_star.stellar_mass)

        disk_codes.append(s_code)
        disk_codes_indices[new_star.key[0]] = len(disk_codes) - 1
        diverged_disks[s_code] = False

        # Saving these values to keep track of dispersed disks later on
        new_star.dispersed_disk_mass = 0.01 * new_star.disk_mass  # Disk is dispersed if it has lost 99% of its initial mass
        new_star.dispersion_threshold = 1E-5 | units.g / units.cm ** 2  # Density threshold for dispersed disks, Ingleby+ 2009
        new_star.dispersed = False
        new_star.checked = False  # I need this to keep track of dispersed disk checks
        new_star.dispersal_time = time
        new_star.photoevap_mass_loss = 0 | units.MJupiter
        new_star.cumulative_photoevap_mass_loss = 0 | units.MJupiter
        new_star.truncation_mass_loss = 0 | units.MJupiter
        new_star.cumulative_truncation_mass_loss = 0 | units.MJupiter
        new_star.EUV = False  # For photoevaporation regime  # TODO I think this is not needed anymore
        new_star.nearby_supernovae = False

        # Copying these back because the disk radius and mass differ slightly from the actual numbers given
        new_star.disk_radius = new_disk.get_radius()
        new_star.disk_mass = new_disk.get_mass(new_star.disk_radius)

        # Initial values of disks
        new_star.initial_disk_size = new_star.disk_radius
        new_star.initial_disk_mass = new_star.disk_mass
        print "disk radius: ", new_star.disk_radius.in_(units.au), new_disk.get_radius()
        print "disk mass: ", new_star.disk_mass.in_(units.MSun), new_disk.get_mass(new_disk.get_radius())

        # Value to keep track of disk sizes and masses as not influenced by photoevaporation
        new_star.disk_size_np = new_star.initial_disk_size
        new_star.disk_mass_np = new_star.initial_disk_mass

    else:
        new_star.disk_radius = 0 | units.au
        new_star.disk_mass = 0 | units.MSun
        new_star.code = False

    return new_star, delay_time


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


def column_density(grid,
                   rc,
                   mass,
                   lower_density=1E-12 | units.g / units.cm**2):
    """ Disk column density definition as in Eqs. 1, 2, and 3 of the paper.
        (Lynden-Bell & Pringle, 1974: Anderson et al. 2013)

    :param grid: disk grid
    :param rc: characteristic disk radius
    :param mass: disk mass
    :param lower_density: density limit for defining disk edge
    :return: disk column density in g / cm**2
    """
    r = grid.value_in(units.au) | units.au
    rd = rc  # Anderson et al. 2013
    Md = mass

    Sigma_0 = Md / (2 * numpy.pi * rc ** 2 * (1 - numpy.exp(-rd / rc)))
    Sigma = Sigma_0 * (rc / r) * numpy.exp(-r / rc) * (r <= rc) + lower_density
    return Sigma


def distance(star1,
             star2):
    """ Return distance between star1 and star2

    :param star1: AMUSE particle
    :param star2: AMUSE particle
    :return: distance in units.parsec
    """
    return numpy.sqrt((star2.x - star1.x)**2 + (star2.y - star1.y)**2 + (star2.z - star1.z)**2)


def luminosity_fit(mass):
    """
    Return stellar luminosity (in LSun) for corresponding mass, as calculated with Martijn's fit

    :param mass: stellar mass in MSun
    :return: stellar luminosity in LSun
    """
    if 0.12 < mass < 0.24:
        return (1.70294E16 * numpy.power(mass, 42.557)) | units.LSun
    elif 0.24 < mass < 0.56:
        return (9.11137E-9 * numpy.power(mass, 3.8845)) | units.LSun
    elif 0.56 < mass < 0.70:
        return (1.10021E-6 * numpy.power(mass, 12.237)) | units.LSun
    elif 0.70 < mass < 0.91:
        return (2.38690E-4 * numpy.power(mass, 27.199)) | units.LSun
    elif 0.91 < mass < 1.37:
        return (1.02477E-4 * numpy.power(mass, 18.465)) | units.LSun
    elif 1.37 < mass < 2.07:
        return (9.66362E-4 * numpy.power(mass, 11.410)) | units.LSun
    elif 2.07 < mass < 3.72:
        return (6.49335E-2 * numpy.power(mass, 5.6147)) | units.LSun
    elif 3.72 < mass < 10.0:
        return (6.99075E-1 * numpy.power(mass, 3.8058)) | units.LSun
    elif 10.0 < mass < 20.2:
        return (9.73664E0 * numpy.power(mass, 2.6620)) | units.LSun
    elif 20.2 < mass:
        return (1.31175E2 * numpy.power(mass, 1.7974)) | units.LSun
    else:
        return 0 | units.LSun


def radiation_at_distance(rad, d):
    """ Return radiation rad at distance d

    :param rad: total radiation from star in erg/s
    :param d: distance in cm
    :return: radiation of star at distance d, in erg * s^-1 * cm^-2
    """
    return rad / (4 * numpy.pi * d**2) | (units.erg / (units.s * units.cm**2))


def FUV_radiation_on_star(star, bright_stars):
    """ Calculate the total FUV radiation over a low mass star, in units of G0

    :param star: star (AMUSE particle) to calculate radiation on
    :param bright_stars: AMUSE particle set of radiating (high mass) stars
    :return: total FUV radiation in G0 units and EUV=True is star if also subject to EUV photoevaporation
    """
    total_radiation = 0
    EUV = False

    for bs in bright_stars:  # For each massive/bright star
        # Calculate FUV luminosity of the bright star, in LSun
        lum = luminosity_fit(bs.stellar_mass.value_in(units.MSun))

        # Calculate distance to bright star
        dist = distance(bs, star)

        # EUV regime -- Use Johnstone, Hollenbach, & Bally 1998
        dmin = 5. * 1E17 * 0.25 * numpy.sqrt(star.disk_radius.value_in(units.cm) / 1E14) | units.cm

        if dist < dmin:
            EUV = True

        else:
            # Other bright stars can still contribute FUV radiation
            radiation_ss = radiation_at_distance(lum.value_in(units.erg / units.s),
                                                 dist.value_in(units.cm)
                                                 )

            radiation_ss_G0 = radiation_ss.value_in(units.erg / (units.s * units.cm ** 2)) / 1.6E-3
            total_radiation += radiation_ss_G0

    return total_radiation, EUV


def photoevaporation_mass_loss(star, bright_stars, dt):
    global FRIED_grid, grid_log10Mdot, grid_stellar_mass, grid_FUV, grid_disk_mass, grid_disk_radius

    radiation_G0, EUV = FUV_radiation_on_star(star, bright_stars)

    if EUV:
        # Photoevaporative mass loss in MSun/yr. Eq 20 from Johnstone, Hollenbach, & Bally 1998
        # From the paper: e ~ 3, x ~ 1.5
        photoevap_Mdot = 2. * 1E-9 * 3 * 4.12 * (star.disk_radius.value_in(units.cm) / 1E14)

        # Calculate total mass lost due to EUV photoevaporation during dt, in MSun
        total_photoevap_mass_loss_euv = float(photoevap_Mdot * dt.value_in(units.yr)) | units.MSun
    else:
        total_photoevap_mass_loss_euv = 0.0 | units.MSun

    # FUV regime -- Use FRIED grid

    # For the small star, I want to interpolate the photoevaporation mass loss
    # xi will be the point used for the interpolation. Adding star values...
    xi = numpy.ndarray(shape=(1, 4), dtype=float)
    xi[0][0] = star.stellar_mass.value_in(units.MSun)
    xi[0][1] = radiation_G0
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
    FUV_i, FUV_j = find_indices(grid_FUV, radiation_G0)
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
                #gravity_sinks = Gravity(ph4, local_sinks)
                #gravity_sinks_to_framework = gravity_sinks.code.particles.new_channel_to(local_sinks)
                #framework_to_gravity_sinks = local_sinks.new_channel_to(gravity_sinks.code.particles)

                converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
                gravity_sinks = ph4(converter, number_of_workers=12)
                gravity_sinks.parameters.timestep_parameter = 0.01
                gravity_sinks.parameters.epsilon_squared = (100 | units.au) ** 2
                gravity_sinks.particles.add_particles(local_sinks)

                gravity_sinks_to_framework = gravity_sinks.particles.new_channel_to(local_sinks)
                framework_to_gravity_sinks = local_sinks.new_channel_to(gravity_sinks.particles)


                #break
                # From this point on I will simply keep evolving the gravity and gravity_sinks codes

        for sink in local_sinks:
            # Iterate over local_sinks instead of hydro.sink_particles so that the leftover
            # sinks can still form stars after the gas code is stopped.

            if sink.mass > IMF_masses[current_mass] and sink.form_star:
                stars_from_sink, delay_time = make_star_from_sink(sink, IMF_masses[current_mass], time)
                sink.form_star = False
                sink.time_threshold = time + delay_time  # Next time at which this sink should form a star

                # Update sink masses and delay times from framework to codes
                if gravity_sinks is None:
                    framework_to_hydro_sinks.copy()
                else:
                    framework_to_gravity_sinks.copy()

                current_mass += 1
                stars.add_particles(stars_from_sink)

                # To keep track of "disked" and "FUV-radiating" stars
                if stars_from_sink.bright:
                    high_mass.append(stars_from_sink.key)
                else:
                    low_mass.append(stars_from_sink.key)  # TODO ojo: Mstar < 0.05 MSun will not have disks


                # I don't care about gravity_sinks for this block
                if gravity is None:
                    print 'Starting gravity code and Bridge'
                    gravity_offset_time = time
                    #gravity = Gravity(ph4, stars, number_of_workers=12)

                    converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
                    gravity = ph4(converter, number_of_workers=12)
                    gravity.parameters.timestep_parameter = 0.01
                    gravity.parameters.epsilon_squared = (100 | units.au) ** 2
                    gravity.particles.add_particles(stars_from_sink)

                    # Enable stopping condition for dynamical encounters
                    dynamical_encounter = gravity.stopping_conditions.collision_detection
                    dynamical_encounter.enable()

                    #gravity_to_framework = gravity.code.particles.new_channel_to(stars)
                    gravity_to_framework = gravity.particles.new_channel_to(stars)
                    #framework_to_gravity = stars.new_channel_to(gravity.code.particles)
                    framework_to_gravity = stars.new_channel_to(gravity.particles)
                    gravity_to_framework.copy()

                    gravhydro = Bridge()
                    gravhydro.add_system(gravity, (hydro,))
                    gravhydro.add_system(hydro, (gravity,))
                    gravhydro.timestep = 0.1 * dt
                else:
                    #gravity.code.particles.add_particles(stars_from_sink)
                    gravity.particles.add_particles(stars_from_sink)
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

            # TODO also here add stellar evolution, disk evolution, photoevap, etc

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
    print stars.disk_radius.in_(units.au)
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

    SFE = 0.01
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

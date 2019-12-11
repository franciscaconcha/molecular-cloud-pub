import numpy

from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from amuse.couple.bridge import Bridge
from amuse.ic.fractalcluster import new_fractal_cluster_model

from cooling_class import SimplifiedThermalModel, SimplifiedThermalModelEvolver
from hydrodynamics_class import Hydro
from gravity_class import Gravity


def write_data(path, timestamp, hydro, index=0, stars=Particles(0)):
    hydro.write_set_to_file(path, index=index)
    filename = "{0}/hydro_stars_particles_i{1:04}.amuse".format(path, index)
    if len(stars) > 0:
        write_set_to_file(stars, filename, "hdf5",
                          timestamp=timestamp,
                          append_to_file=False)


def fill_mass_function_with_sink_mass(total_mass, method):
    # print "Make mass function for M=", total_mass.in_(units.MSun)
    masses = [] | units.MSun

    while total_mass > 0 | units.MSun:
        mass = new_kroupa_mass_distribution(1, mass_max=100 | units.MSun)[0]
        if mass > total_mass:
            mass = total_mass
        total_mass -= mass
        masses.append(mass)

    return masses


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
    stars_from_sink.mass = stellar_mass
    sink.mass -= stars_from_sink.mass

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
                gravity_sinks.evolve_model(time - gravity_offset_time)
                print gravity_sinks.particles
                gravity_to_framework.copy()
                gravity_sinks_to_framework.copy()
                #local_sinks.synchronize_to(hydro.sink_particles)

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
            write_data(save_path, time, hydro=hydro, index=index, stars=stars)

            # Saving stars and local sinks separately from the Hydro files
            # Saving star particles
            write_set_to_file(stars,
                              '{0}/stars_{1}Myr.hdf5'.format(save_path,
                                                             time.value_in(units.Myr)),
                              'hdf5')
            # Saving local sink particles
            write_set_to_file(local_sinks,
                              '{0}/sinks_{1}Myr.hdf5'.format(save_path,
                                                             time.value_in(units.Myr)),
                              'hdf5')

    #print len(gravity.code.particles)
    hydro.stop()
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

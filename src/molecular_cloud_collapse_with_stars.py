import numpy

from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from cooling_class import SimplifiedThermalModel, SimplifiedThermalModelEvolver
from hydrodynamics_class import Hydro
from amuse.couple.bridge import Bridge


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


def make_star_from_sink(sink, stellar_mass, time, factor):

    # Delay time for next star formation is a decay from the tff of the sink
    delay_time = sink.time_threshold * numpy.exp(-time.value_in(units.Myr))

    print "Forming star of mass {0} from sink mass {1}".format(stellar_mass.in_(units.MSun),
                                                               sink.mass.in_(units.MSun))

    new_star = Particles(1)
    new_star.stellar_mass = stellar_mass
    new_star.mass = stellar_mass

    if new_star.stellar_mass <= 1.9 | units.MSun:
        sink.mass -= 1.1 * new_star.stellar_mass  # Also removing mass of the future disk
    else:
        sink.mass -= new_star.stellar_mass  # No disk here

    # 'Normal' star parameters: location, velocity, etc
    # Find position offset inside sink radius
    Rsink = sink.radius.value_in(units.parsec)
    offsetx = numpy.random.uniform(-factor * Rsink, factor * Rsink) | units.parsec
    offsety = numpy.random.uniform(-factor * Rsink, factor * Rsink) | units.parsec
    offsetz = numpy.random.uniform(-factor * Rsink, factor * Rsink) | units.parsec
    print offsetx, offsety, offsetz
    new_star.x = sink.x + offsetx
    new_star.y = sink.y + offsety
    new_star.z = sink.z + offsetz

    new_star.vx = sink.vx
    new_star.vy = sink.vy
    new_star.vz = sink.vz

    new_star.tborn = time

    return new_star, delay_time


def run_molecular_cloud(gas_particles, sink_particles, SFE, tstart, tend, dt_diag, save_path, factor, index=0):

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
    IMF_masses = new_kroupa_mass_distribution(10000,
                                              mass_max=150 | units.MSun,
                                              random=True)  # Randomized order
    print "Some masses, pre >0.08:"
    print IMF_masses[0], IMF_masses[100], IMF_masses[200]
    IMF_masses = [m for m in IMF_masses if m >= 0.08 | units.MSun]  # Wall+2019 stellar mass range
    current_mass = 0  # To keep track of formed stars

    print "Some masses, post >0.08:"
    print IMF_masses[0], IMF_masses[100], IMF_masses[200]

    stars = Particles(0)  # Here we keep the newly formed stars

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
            if len(hydro.sink_particles) > 0:
                Mtot = hydro.gas_particles.mass.sum() + hydro.sink_particles.mass.sum()
            else:
                Mtot = hydro.gas_particles.mass.sum()

            if gravity is not None:
                Mtot += gravity.particles.mass.sum()

            if len(local_sinks) > 0:
                MC_SFE = local_sinks.mass.sum() / Mtot
            else:
                MC_SFE = 0

            if MC_SFE >= SFE:
                print "************** SFE reached, sinks will stop forming **************"
                sink_formation = False

                # Final synchro and copy of hydro sinks to local sinks

                hydro.sink_particles.synchronize_to(local_sinks)
                hydro_sinks_to_framework.copy()

                converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)

                # This gravity_sinks code was to keep the sinks in the gravity code, so that they
                # gravitationally interact with the stars. Won't be using it anymore (for now...)
                gravity_sinks = ph4(converter)#, number_of_workers=12)
                gravity_sinks.parameters.timestep_parameter = 0.01
                gravity_sinks.parameters.epsilon_squared = (100 | units.au) ** 2
                gravity_sinks.particles.add_particles(local_sinks)

                gravity_sinks_to_framework = gravity_sinks.particles.new_channel_to(local_sinks)
                framework_to_gravity_sinks = local_sinks.new_channel_to(gravity_sinks.particles)

                framework_to_gravity_sinks.copy()

                gravity_time_offset = time

            else:
                if gravhydro is None:
                    hydro.evolve_model(time)
                else:
                    print "EVOLVING GRAVHYDRO"
                    gravhydro.evolve_model(time)
                    gravity_to_framework.copy()

                if len(local_sinks):
                    print "before synch:"
                    print local_sinks.mass.in_(units.MSun)
                # Synchronize sinks, then update local sinks
                hydro.sink_particles.synchronize_to(local_sinks)
                hydro_sinks_to_framework.copy()

                if len(local_sinks):
                    print "after synch:"
                    print local_sinks.mass.in_(units.MSun)

        else:  # sink_formation == False
            if gravity_sinks is None:
                #print "No gravity_sinks yet"
                pass
            else:
                gravity_sinks.evolve_model(time - gravity_time_offset)
                gravity_sinks_to_framework.copy()

            if gravity is None:
                hydro.evolve_model(time)
                #print "No gravity yet"
            else:
                print "EVOLVING GRAVITY ONLY"
                gravity.evolve_model(time)
                gravity_to_framework.copy()
                # no need to update local_sinks now because i dont want them to keep accreting

        # This loop will be executed regardless if the hydro code is active or not (sink_formation)
        print "Trying to form a star of mass ", IMF_masses[current_mass]
        print "N hydro sinks = {0}".format(len(hydro.sink_particles))
        print "N local_sinks = {0}".format(len(local_sinks))
        for sink in local_sinks:
            # Iterate over local_sinks instead of hydro.sink_particles so that the leftover
            # sinks can still form stars after the gas code is stopped.

            if sink.mass >= IMF_masses[current_mass] and sink.form_star:
                # Make a star!
                print "Making a star"
                #print "Sink mass before: ", sink.mass.value_in(units.MSun)
                stars_from_sink, delay_time = make_star_from_sink(sink, IMF_masses[current_mass], time, factor)
                sink.form_star = False
                sink.time_threshold = time + delay_time  # Next time at which this sink should form a star

                # Update sink masses and delay times from framework to codes
                if gravity_sinks is None:
                    framework_to_hydro_sinks.copy()
                else:
                    framework_to_gravity_sinks.copy()

                if gravity is None:
                    converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
                    gravity_offset_time = time
                    gravity = ph4(converter)#, number_of_workers=12)
                    gravity.particles.add_particles(stars_from_sink)
                    gravity_to_framework = gravity.particles.new_channel_to(stars)
                    framework_to_gravity = stars.new_channel_to(gravity.particles)
                    gravhydro = Bridge()
                    gravhydro.add_system(gravity, (hydro,))
                    gravhydro.add_system(hydro, (gravity,))
                    gravhydro.timestep = 0.1 * dt
                else:
                    gravity.particles.add_particles(stars_from_sink)
                    gravity_to_framework.copy()

                current_mass += 1
                stars.add_particles(stars_from_sink)
                print "stars = {0}, gravity.particles = {1}".format(len(stars),
                                                                    len(gravity.particles))
                #print "Sink mass after: ", sink.mass.value_in(units.MSun)

            elif sink.mass >= IMF_masses[current_mass] and not sink.form_star:
                print "Sink is massive enough, but it's not yet time to form a star."
                if time >= sink.time_threshold:
                    sink.form_star = True

                    # Update sink.form_star in corresponding code
                    #if gravity_sinks is None:
                    framework_to_hydro_sinks.copy()
                    #else:
                    #    framework_to_gravity_sinks.copy()
                    # sink will form a star in the next timestep

            elif sink.mass < IMF_masses[current_mass] and sink.form_star:
                print "Sink is not massive enough to form this star."
                sinks_masses = [s.mass for s in local_sinks]

                # No sink has enough mass to form this star
                if not sink_formation and all(i < IMF_masses[current_mass] for i in sinks_masses):
                    # If the sinks are not accreting anymore, I keep moving through the IMF
                    # so that we don't get stuck trying to form a massive star with low mass sinks
                    current_mass += 1
                # sink.form_star = False

            # Stopping condition for star formation -- if this is reached the code finishes
            sinks_masses = [s.mass for s in local_sinks]
            if len(IMF_masses[current_mass:]) > 0:
                if not sink_formation and all(i < min(IMF_masses[current_mass:]) for i in sinks_masses):
                    print "All mass in sinks has ran out -- all possible stars have been formed!"
                    print "Made {0} stars.".format(len(stars))
                    print tend, time
                    time += tend  # To break out of outer loop
                    break
            else:
                print "All stars in the IMF have been formed."
                print "Made {0} stars.".format(len(stars))
                print "Finishing at {0}".format(time.in_(units.Myr))
                print tend, time
                time += tend  # To break out of outer loop
                break
            #gravity_to_framework.copy()
            #if gravity_sinks is None:
            framework_to_hydro_sinks.copy()
            #else:
            #    framework_to_gravity_sinks.copy()

        E = hydro.gas_particles.kinetic_energy() \
            + hydro.gas_particles.potential_energy() \
            + hydro.gas_particles.thermal_energy()
        E_th = hydro.gas_particles.thermal_energy()
        Eerr = (E - E0) / E0
        print 'energy=', E, 'energy_error=', Eerr, 'e_th=', E_th
        print "maximal_density:", gas_particles.rho.max().in_(units.MSun / units.parsec ** 3)
        #for s in local_sinks:
        #    print s.mass.value_in(units.MSun)

        #hydro.print_diagnostics()
        #if gravhydro is None:
        #    print "No gravhydro yet."
        #else:
        #    print "gravhydro"
            # print_diagnostics(gravhydro)
        if time > t_diag:
            index += 1
            t_diag += dt

            if sink_formation:
                write_data(save_path, time, hydro=hydro, index=index, stars=stars)
            else:
                gravity_sinks_to_framework.copy()
                gravity_to_framework.copy()
                write_data(save_path, time, hydro=hydro, index=index, stars=stars)
                # Saving local sink particles
                write_set_to_file(local_sinks,
                                  '{0}/hydro_sink_particles_i00{1}.amuse'.format(save_path,
                                                                                 index),
                                  # Because I add tend to break out of loop!
                                  'hdf5')

    # Saving stars and local sinks separately from the Hydro files
    # Saving star particles
    write_set_to_file(stars,
                      '{0}/gravity_stars_t{1:.2f}Myr.hdf5'.format(save_path,
                                                                  (time - tend).value_in(units.Myr)), # Because I add tend to break out of loop!
                      'hdf5')
    # Saving local sink particles
    write_set_to_file(local_sinks,
                      '{0}/gravity_sinks_t{1:.2f}Myr.hdf5'.format(save_path,
                                                                  (time - tend).value_in(units.Myr)), # Because I add tend to break out of loop!
                      'hdf5')

    #print len(gravity.code.particles)
    hydro.stop()#, gravity.stop(), gravity_sinks.stop()
    print "Made {0} stars.".format(len(stars))
    #print stars.disk_radius.in_(units.au)
    return gas_particles


def main(filename, save_path, tend, dt_diag, Ncloud, Mcloud, Rcloud, factor):
    import datetime
    print("START: {0}".format(datetime.datetime.now()))

    if len(filename) == 0:  # Start a new run
        try:
            import os
            os.makedirs(save_path)
        except OSError, e:
            if e.errno != 17:
                raise
            pass

        gas_particles = generate_initial_conditions_for_molecular_cloud(Ncloud,
                                                                        Mcloud=Mcloud,
                                                                        Rcloud=Rcloud)
        rho_cloud = Mcloud / Rcloud ** 3
        tff = 0.5427 / numpy.sqrt(constants.G * rho_cloud)
        print "Freefall timescale=", tff.in_(units.Myr)
        rho_cloud = 3. * Mcloud / (4. * numpy.pi * Rcloud ** 3)
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
        print "Starting with {0} sinks".format(len(sink_particles))
    else:
        sink_particles = Particles(0)

    print "Time= {0}".format(start_time.in_(units.Myr))
    print "index = {0}, Ngas = {1}, Nsinks = {2}".format(index, len(gas_particles), len(sink_particles))

    SFE = 0.3  # Star formation efficiency 30%

    parts = run_molecular_cloud(gas_particles, sink_particles, SFE, start_time, tend, dt_diag, save_path, factor, index)

    print("END: {0}".format(datetime.datetime.now()))


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
    result.add_option("--fcloud", dest="factor",
                      type="int",
                      default=1,
                      help="factor for Rsink for location of new stars")

    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)

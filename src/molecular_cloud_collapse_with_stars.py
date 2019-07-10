import numpy
from matplotlib import pyplot 

from amuse.lab import *
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from amuse.couple.bridge import Bridge


from cooling_class import SimplifiedThermalModel, SimplifiedThermalModelEvolver
from hydrodynamics_class import Hydro
from gravity_class import Gravity

def write_data(path, hydro, index=0, stars=Particles(0)):
        hydro.write_set_to_file(path, index=index)
        filename = "{0}/hydro_stars_particles_i{1:04}.amuse".format(path, index)
        if len(stars)>0:
                write_set_to_file(stars, filename, "hdf5",
                                  timestamp=hydro.model_time,
                                  append_to_file=False)

def fill_mass_function_with_sink_mass(total_mass):
    print "Make mass function for M=", total_mass.in_(units.MSun)
    masses = new_kroupa_mass_distribution(1000, 100 | units.MSun)
    mean_mass = masses.sum() / len(masses)
    Nstars = int(total_mass/mean_mass)
    masses = new_kroupa_mass_distribution(Nstars, 100|units.MSun)
    print "Generate N=", Nstars, "stars of mean mass m=", mean_mass.in_(units.MSun)
    print "total mass generated M=", masses.sum().in_(units.MSun)
    return masses


def print_diagnostics(gravhydro):
    print "Time =", gravhydro.model_time.in_(units.Myr)
    print "Ngravhydro =", len(gravhydro.particles), "Mgravhydro=", gravhydro.particles.mass.sum().in_(units.MSun)
    print " "


def generate_initial_conditions_for_molecular_cloud(N, Mcloud, Rcloud):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)
    gas = molecular_cloud(targetN=N,convert_nbody=conv,
            base_grid=body_centered_grid_unit_cube).result
    gas.name = "gas"

    return gas


def run_molecular_cloud(gas_particles, sink_particles, tstart, tend, dt_diag, save_path, index=0):

    Mcloud = gas_particles.mass.sum()
    
    stars = Particles(0)
    hydro = Hydro(Fi, gas_particles, tstart)

    if len(sink_particles) > 0:
        hydro.sink_particles.add_particles(sink_particles)
        hydro.sink_particles.synchronize_to(hydro.code.dm_particles)


    gravity = None
    gravhydro = None
    dt = min(dt_diag, 0.1 | units.Myr)
    t_diag = 0 | units.Myr

    mass_treshold_for_star_formation = 10 | units.MSun    
    E0 = hydro.gas_particles.kinetic_energy() \
         + hydro.gas_particles.potential_energy() \
         + hydro.gas_particles.thermal_energy()
    time = gas_particles.get_timestamp()

    while time < tend:
        time += dt
        print "Evolve to time=", time.in_(units.Myr)
        Mtot = 0|units.MSun
                
        if len(hydro.sink_particles) > 0:
            print "Mass conservation at t = {0}:".format(time.in_(units.Myr))
            print "Slocal_gas = {0}, " \
                  "Slocal_sinks = {1}, " \
                  "sum = {2}".format(hydro.gas_particles.mass.sum().in_(units.MSun),
                                      hydro.sink_particles.mass.sum().in_(units.MSun),
                                      (hydro.gas_particles.mass.sum() + hydro.sink_particles.mass.sum()).in_(units.MSun))
            print "Shydro_gas = {0}, " \
                  "Shydro_sinks = {1}, " \
                  "sum = {2}\n*".format(hydro.code.gas_particles.mass.sum().in_(units.MSun),
                                        hydro.code.dm_particles.mass.sum().in_(units.MSun),
                                        (hydro.code.gas_particles.mass.sum() + hydro.code.dm_particles.mass.sum()).in_(units.MSun))

            Mtot = hydro.gas_particles.mass.sum() + hydro.sink_particles.mass.sum()

            removed_sinks = Particles(0)

            for sink in hydro.sink_particles:
                if sink.mass > mass_treshold_for_star_formation:
                    print "Turn sink into cluster. Msink = {0}".format(sink.mass.in_(units.MSun))
                    masses = fill_mass_function_with_sink_mass(sink.mass)
                    local_converter = nbody_system.nbody_to_si(masses.sum(), sink.radius)
                    stars_from_sink = new_plummer_model(len(masses), local_converter)
                    stars_from_sink.mass = masses
                    stars_from_sink.scale_to_standard(local_converter)
                    stars_from_sink.age = time
                    removed_sinks.add_particle(sink)

                    stars.add_particles(stars_from_sink)
                    if gravity==None:
                        gravity_offset_time = time
                        gravity = Gravity(ph4, stars)
                        gravhydro = Bridge()
                        gravhydro.add_system(gravity, (hydro.code,))
                        gravhydro.add_system(hydro.code, (gravity,))
                        gravhydro.timestep = 0.1*dt
                    else:
                        gravity.code.particles.add_particles(stars_from_sink)
                        gravity.update_channels(stars)

            if len(removed_sinks)>0:
                # clean up hydro code by removing sink particles.
                print "Clean up hydro code by removing sink particles."
                hydro.sink_particles.remove_particle(removed_sinks)
                hydro.sink_particles.synchronize_to(hydro.code.dm_particles)
                    
        else:
            print "Mass conservation at t = {0}:".format(time.in_(units.Myr))
            print "Local: {0}, Hydro: {1}\n".format(hydro.gas_particles.mass.sum().in_(units.MSun),
                                                            hydro.code.gas_particles.mass.sum().in_(units.MSun))
            Mtot = hydro.gas_particles.mass.sum()

        if Mtot < Mcloud - (1.e-5|units.MSun):
            print "Mass is not conserved: Mtot = {0} MSun, Mcloud = {1} MSun".format(Mtot.in_(units.MSun),
                                                                                     Mcloud.in_(units.MSun))
            exit(-1)

        if gravhydro == None:
            hydro.evolve_model(time)
        else:
            gravhydro.evolve_model(time)
            
        E = hydro.gas_particles.kinetic_energy() \
             + hydro.gas_particles.potential_energy() \
             + hydro.gas_particles.thermal_energy()
        E_th = hydro.gas_particles.thermal_energy()
        Eerr = (E-E0)/E0
        print 'energy=', E, 'energy_error=', Eerr, 'e_th=', E_th
        print "maximal_density:",gas_particles.rho.max().in_(units.MSun/units.parsec**3)

        hydro.print_diagnostics()
        if gravhydro==None:
            print "No gravhydro yet."
        else:
            print_diagnostics(gravhydro)
        if time>t_diag:
            index=index+1
            t_diag += dt_diag
            write_data(save_path, hydro=hydro, index=index, stars=stars)

    hydro.stop()
    return gas_particles


def main(filename, save_path, tend, dt_diag, Ncloud, Mcloud, Rcloud):

    if len(filename) == 0:
        gas_particles = generate_initial_conditions_for_molecular_cloud(o.Ncloud,
                                                                        Mcloud=o.Mcloud,
                                                                        Rcloud=o.Rcloud)
        rho_cloud = o.Mcloud / o.Rcloud ** 3
        tff = 0.5427 / numpy.sqrt(constants.G * rho_cloud)
        print "Freefall timescale=", tff.in_(units.Myr)
        rho_cloud = 3. * o.Mcloud / (4. * numpy.pi * o.Rcloud ** 3)
        print rho_cloud
        tend = 10 * tff
        dt_diag = 0.1 * tff
        hydro = Hydro(Fi, gas_particles)
        write_data(save_path, hydro)
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

    parts = run_molecular_cloud(gas_particles, sink_particles, start_time, tend, dt_diag, save_path, index)

  
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename",
                      default="",
                      help="input filename")
    result.add_option("-s", dest="save_path",
                      default=".",
                      help="save path for results")
    result.add_option("--tend", dest="tend",
                      unit=units.Myr,
                      type="float",
                      default = 1.0|units.Myr,
                      help="end time")
    result.add_option("--dt_diag", dest="dt_diag",
                      unit=units.Myr,
                      type="float",
                      default=0.1|units.Myr,
                      help="diagnosticstime step")
    result.add_option("--Ncloud", dest="Ncloud",
                      default=1000,
                      type="float",
                      help="number of gas particles.")
    result.add_option("--Mcloud", dest="Mcloud",
                      unit=units.MSun,
                      type="float",
                      default=1000|units.MSun,
                      help="cloud mass")
    result.add_option("--Rcloud", dest="Rcloud",
                      unit=units.parsec,
                      type="float",
                      default=3|units.parsec,
                      help="cloud size")

    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    numpy.random.seed(3141)
    main(**o.__dict__)

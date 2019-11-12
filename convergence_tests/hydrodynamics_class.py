import time
import numpy

from amuse.lab import *
from amuse.units import units, nbody_system
from amuse.datamodel import Particles

from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2

from cooling_class import Cooling, SimplifiedThermalModelEvolver
# from amuse.ext.sink import SinkParticles
from amuse.ext.sink import new_sink_particles

COOL = True


class Hydro:
    sink_id = 0
    g_newsinks = 0

    def __init__(self, hydro_code, particles, tstart=0 | units.Myr):

        if not hydro_code in [Fi, Gadget2]:
            raise Exception("unsupported Hydro code: %s" % (hydro_code.__name__))

        self.current_time = tstart
        self.typestr = "Hydro"
        self.namestr = hydro_code.__name__

        system_size = 2 | units.parsec

        eps = 0.05 | units.parsec
        Mcloud = particles.mass.sum()
        N = len(particles)
        dt = 0.004 * numpy.pi * numpy.power(eps, 1.5) / numpy.sqrt(constants.G * Mcloud / N)
        print "Hydro timesteps:", dt, "N=", len(particles)

        self.gas_particles = particles
        self.sink_particles = Particles(0)

        self.cooling_flag = "thermal_model"

        self.density_threshold = (1 | units.MSun) / (eps) ** 3
        self.merge_radius = 0.2 * eps

        self.converter = nbody_system.nbody_to_si(1 | units.MSun, system_size)
        self.sink_attributes = ['name', 'birth_age', 'angular_momentum', 'mass', 'radius', 'x', 'y', 'z', 'vx', 'vy',
                                'vz', 'Lx', 'Ly', 'Lz']

        if hydro_code is Fi:
            self.code = hydro_code(self.converter, mode="openmp", redirection="file")

            self.code.parameters.begin_time = 0.0 | units.Myr
            self.code.parameters.use_hydro_flag = True
            self.code.parameters.self_gravity_flag = True
            self.code.parameters.periodic_box_size = 100. * system_size
            """
            isothermal_flag:
            When True then we have to do our own adiabatic cooling (and Gamma has to be 1.0)
            When False then we dont do the adjabatic cooling and Fi is changing u
            """
            self.code.parameters.isothermal_flag = True
            self.code.parameters.integrate_entropy_flag = False
            self.code.parameters.timestep = 0.5 * dt
            # self.code.parameters.verbosity=99
            self.code.parameters.verbosity = 0

            self.code.parameters.gamma = 1
            self.code.parameters.integrate_entropy_flag = False
            self.gamma = self.code.parameters.gamma

        if hydro_code is Gadget2:
            print "WARNING: Gadget support is WIP"
            print "check the Gadget2 makefile_options"
            # todo: variable number_of_workers
            self.code = hydro_code(self.converter, number_of_workers=4)
            self.code.parameters.begin_time = 0.0 | units.Myr
            self.code.parameters.time_max = dt * 2 ** int(numpy.log2(4 * (10 | units.Myr) / dt))
            self.code.parameters.max_size_timestep = 0.5 * dt
            # ~ self.code.parameters.min_size_timestep=
            print "Isofag;", self.code.parameters.isothermal_flag
            print "gamma=", self.code.parameters.polytropic_index_gamma
            # assert self.code.parameters.isothermal_flag == True
            assert self.code.parameters.no_gravity_flag == False
            # constant gas smoothing
            self.code.parameters.gas_epsilon = eps
            assert self.code.parameters.eps_is_h_flag == False
            # assert self.code.parameters.polytropic_index_gamma == 1.

            if self.cooling_flag == "internal":
                raise Exception("gadget internal cooling not implemented")

            self.gamma = self.code.parameters.polytropic_index_gamma

        self.code.parameters.stopping_condition_maximum_density = self.density_threshold
        self.code.commit_parameters()

        if len(self.gas_particles) > 0:
            self.code.gas_particles.add_particles(self.gas_particles)
        if len(self.sink_particles) > 0:
            self.code.dm_particles.add_particles(self.sink_particles)

        self.parameters = self.code.parameters

        print self.code.parameters

        # In order to be using the bridge
        self.get_gravity_at_point = self.code.get_gravity_at_point
        self.get_potential_at_point = self.code.get_potential_at_point
        self.get_hydro_state_at_point = self.code.get_hydro_state_at_point

        # Create a channel
        self.channel_to_gas = self.code.gas_particles.new_channel_to(self.gas_particles)
        self.channel_to_sinks = self.code.dm_particles.new_channel_to(self.sink_particles)
        self.channel_from_sinks = self.sink_particles.new_channel_to(self.code.dm_particles)

        # External Cooling
        print "Cooling flag:", self.cooling_flag
        self.cooling = SimplifiedThermalModelEvolver(self.code.gas_particles)
        # self.cooling = Cooling(self.code.gas_particles)
        self.cooling.model_time = self.model_time

    def print_diagnostics(self):
        print "Time=", self.model_time.in_(units.Myr)
        print "N=", len(self.gas_particles), len(self.sink_particles)

        if len(self.sink_particles) > 0:
            print "Sink masses:", len(self.code.dm_particles.mass)
            print "Sink masses:", len(self.sink_particles.mass)

    def write_set_to_file(self, save_path, index):
        filename = "{0}/hydro_gas_particles_i{1:04}.amuse".format(save_path, index)
        write_set_to_file(self.gas_particles, filename, "amuse", timestamp=self.model_time,
                          append_to_file=False)
        if len(self.sink_particles):
            filename = "{0}/hydro_sink_particles_i{1:04}.amuse".format(save_path, index)
            write_set_to_file(self.sink_particles, filename, "amuse",
                              timestamp=self.model_time)

    @property
    def model_time(self):
        return self.current_time

    @property
    def gas_particles(self):
        return self.code.gas_particles

    @property
    def stop(self):
        return self.code.stop

    def evolve_model(self, model_time):

        # print "timing:", self.model_time.in_(units.Myr), self.code.model_time.in_(units.Myr), self.current_time.in_(units.Myr)

        start_time = time.time()

        density_limit_detection = self.code.stopping_conditions.density_limit_detection
        density_limit_detection.enable()

        # this needs to be cleaned up.
        model_time_old = self.model_time
        self.current_time = model_time
        dt = model_time - model_time_old
        hydro_time = self.code.model_time + dt
        print "Evolve Hydrodynamics:", dt.in_(units.Myr), self.model_time.in_(units.Myr), self.code.model_time.in_(
            units.Myr), self.current_time.in_(units.Myr), hydro_time.in_(units.Myr)

        print "Density treshold:", self.density_threshold.in_(
            units.g / units.cm ** 3), self.code.gas_particles.density.max().in_(units.g / units.cm ** 3)

        if COOL:
            print "Cool gas for dt=", (dt / 2).in_(units.Myr)
            self.cooling.evolve_for(dt / 2)
            # print "...done."
        self.code.evolve_model(hydro_time)

        # print "gas evolved."
        while density_limit_detection.is_set():
            self.resolve_sinks()

            print "..done"
            self.code.evolve_model(hydro_time)
            self.channel_to_sinks.copy()
            print "end N=", len(self.sink_particles), len(self.code.dm_particles)

        if COOL:
            print "Cool gas for another dt=", (dt / 2).in_(units.Myr)
            self.cooling.evolve_for(dt / 2)
            # print "...done."

        self.merge_sinks()

        if len(self.sink_particles) > 0:
            sinks = new_sink_particles(self.sink_particles)
            sinks.accrete(self.gas_particles)
            for si in range(len(self.sink_particles)):
                self.sink_particles[si].Lx += sinks[si].angular_momentum[0]
                self.sink_particles[si].Ly += sinks[si].angular_momentum[1]
                self.sink_particles[si].Lz += sinks[si].angular_momentum[2]

            self.gas_particles.synchronize_to(self.code.gas_particles)
            # make sure that the accreted mass is copied to the Hydro code..+++
            self.channel_from_sinks.copy()

        # self.code.evolve_model(model_time)
        print "final N=", len(self.sink_particles), len(self.code.dm_particles)
        print "final Ngas=", len(self.gas_particles), len(self.code.gas_particles)
        self.channel_to_gas.copy()
        self.channel_to_sinks.copy()
        print "final N=", len(self.sink_particles), len(self.code.dm_particles)
        if len(self.sink_particles) > 0:
            print "mean sink mass:", self.sink_particles.mass.min().in_(
                units.MSun), self.sink_particles.mass.mean().in_(units.MSun), self.sink_particles.mass.max().in_(
                units.MSun)
            # ~ print "Hydro arrived at:", self.code.model_time.in_(units.Myr)

            # print "Accrete from ambied gas"
            # self.accrete_sinks_from_ambiant_gas()

    def resolve_sinks(self):
        print "processing high dens particles...",
        highdens = self.gas_particles.select_array(lambda rho: rho > self.density_threshold, ["rho"])
        print "N=", len(highdens)
        print "sinks mass: ", highdens.mass.value_in(units.MSun)
        candidate_sinks = highdens.copy()
        self.gas_particles.remove_particles(highdens)
        self.gas_particles.synchronize_to(self.code.gas_particles)

        print "new sinks..."
        if len(candidate_sinks) > 0:  # had to make some changes to prevent double adding particles
            print "Adding sinks, N=", len(candidate_sinks)
            newsinks_in_code = self.code.dm_particles.add_particles(candidate_sinks)
            newsinks = Particles()
            for nsi in newsinks_in_code:
                if nsi not in self.sink_particles:
                    newsinks.add_particle(nsi)
                else:
                    print "this sink should not exist"
            newsinks.name = "Sink"
            newsinks.birth_age = self.model_time
            newsinks.Lx = 0 | (units.g * units.m ** 2) / units.s
            newsinks.Ly = 0 | (units.g * units.m ** 2) / units.s
            newsinks.Lz = 0 | (units.g * units.m ** 2) / units.s

            for ns in newsinks:
                ns.id = self.sink_id
                ns.merged_ids = numpy.zeros((1, ))
                self.sink_id += 1
                ns.form_star = True
                ns.time_threshold = self.model_time

            self.g_newsinks = len(newsinks)

            print "pre N=", len(self.sink_particles), len(newsinks), len(self.code.dm_particles)
            # self.sink_particles.add_sinks(newsinks)
            self.sink_particles.add_particles(newsinks)
            print "post N=", len(self.sink_particles), len(newsinks), len(self.code.dm_particles)
        else:
            print "N candidates:", len(candidate_sinks)

    def merge_sinks(self):
        if len(self.sink_particles) <= 0:
            return
        print "Let gravity take care of merging sinks"
        if self.sink_particles.radius.max() <= (0 | units.AU):
            return
        print "identify groups.."
        ccs = self.sink_particles.copy().connected_components(threshold=self.merge_radius)
        if len(ccs):
            print "merging sink sets... "
        nmerge = 0
        newsinks = Particles()
        for cc in ccs:
            if len(cc) > 1:
                nmerge += 1
            print "Merge sinks: N= ", len(cc)
            merge_two_sinks(self.sink_particles, cc.copy(), self.sink_id, self.model_time)
            self.sink_id += 2  # To account for the merged sink formed with a new id
            self.sink_particles.synchronize_to(self.code.dm_particles)
            print "sinks merged"


def merge_two_sinks(bodies, particles_in_encounter, id, time):
    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()
    new_particle = Particles(1)
    mu = particles_in_encounter[0].mass / particles_in_encounter.mass.sum()
    new_particle.birth_age = particles_in_encounter.birth_age.min()
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.name = "Sink"
    new_particle.radius = particles_in_encounter.radius.max()

    new_particle.id = id + 1
    new_particle.form_star = True
    new_particle.time_threshold = time

    #l = []
    #for p in particles_in_encounter:
    #    l.append(p.id)
    #new_particle.merged_ids = numpy.pad(numpy.array(l), (0, newsinks - len(particles_in_encounter)), 'constant')



    print "old radius:", particles_in_encounter.radius.value_in(units.AU)
    print "new radius:", new_particle.radius.value_in(units.AU)
    bodies.add_particles(new_particle)
    print "Two sinks (M=", particles_in_encounter.mass, particles_in_encounter.id, ") collided at d=", com_pos.length()
    bodies.remove_particles(particles_in_encounter)

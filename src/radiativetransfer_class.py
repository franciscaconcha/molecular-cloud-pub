import numpy as np

import amuse.lab as al


class Radiative:

    def __init__ (self, rad_code, particles, tstart=0|al.units.Myr):

        if not rad_code in [al.SimpleX]:
            raise Exception('Unsupported Radiative code: %s' % rad_code.__name__)

        self.current_time = tstart
        self.typestr = "Radiative"
        self.namestr = rad_code.__name__

        box_size = 9 | al.units.parsec

        self.particles = particles


        if rad_code is al.SimpleX:

            self.code = rad_code(redirection='file')

            self.code.parameters.blackbody_spectrum_flag = True
            self.code.parameters.source_effective_T = 4.*10.**4. | al.units.K   # ~T of massive stars, but check if that is appropriate
            self.code.parameters.thermal_evolution_flag = True
            self.code.parameters.metal_cooling_flag = True
            self.code.parameters.recombination_radiation_flag = True
            self.code.parameters.collisional_ionization_flag = True
            self.code.parameters.box_size = box_size
            self.code.parameters.hilbert_order = 1


        self.channel_to_particles = self.code.particles.new_channel_to(self.particles)
        self.channel_from_particles = self.particles.new_channel_to(self.code.particles)


    @property
    def model_time (self):
        return self.current_time


    @property
    def particles (self):
        return self.code.particles


    @property
    def stop (self):
        return self.code.stop


    def evolve_model (self, end_time):

        self.channel_to_particles.copy()

        self.code.evolve_model(end_time)

        self.channel_from_particles.copy()

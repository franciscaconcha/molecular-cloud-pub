import os.path
import math
import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.evrard_test import uniform_unit_sphere
import time


from amuse.community.mercury.interface import Mercury
from amuse.community.mikkola.interface import Mikkola
from amuse.units.optparse import OptionParser

from amuse.community.adaptb.interface import Adaptb

class Gravity:
    def __init__(self, gravity_code, stars, converter=None):
        if converter==None:
            self.converter=nbody_system.nbody_to_si(1|units.MSun,1|units.parsec)
        else:
            self.converter = converter
        self.code = gravity_code(self.converter)
        self.code.initialize_code()
        self.get_gravity_at_point = self.code.get_gravity_at_point
        self.get_potential_at_point = self.code.get_potential_at_point
        if hasattr(self.code.parameters, "bs_tolerance"):
            self.code.parameters.bs_tolerance  = 1.0e-10
            self.code.parameters.word_length = 512
        if hasattr(self.code.parameters, "lightspeed"):
            self.code.parameters.lightspeed = 0 | units.kms

        self.code.particles.add_particles(stars)
        if isinstance(self.code, Mercury):
            self.code.commit_particles() 

        self.update_channels(stars)

    @property
    def model_time(self):
        return self.code.model_time
    @property
    def stop(self):
        return self.code.stop
    @property
    def particles(self):
        return self.code.particles

    def update_channels(self, stars):
        self.channel_to_framework = self.code.particles.new_channel_to(stars)
        self.channel_from_framework = stars.new_channel_to(self.code.particles)
        
    def evolve_model(self, model_time):
        for gi in self.channel_from_framework:
            gi.copy()
            #gi.copy_attributes("mass")
        self.code.evolve_model(model_time)
        for gi in self.channel_to_framework:
            gi.copy()

def print_diagnostics(model_time, particles, converter):
    print "Diaganostics: Time=", time, "N=", len(particles)
    
def main(t_end=1687|units.yr, n_steps=1, filename=None):
    black_hole, stars = initialize_sstars(2012|units.yr, S_name, S_a_arcsec, S_ecc, S_inc, S_omra, S_Omega, S_tperi, S_Period)

    gravity = Gravity(Mercury, [black_hole, stars])
#    gravity = Gravity(Mikkola, [black_hole, stars])
    
    print_diagnostics(gravity.model_time, [black_hole, stars], gravity.converter) 
    if filename:
        write_set_to_file(gravity.particles, filename, "hdf5")

    dt = t_end/float(n_steps)
    while gravity.model_time < t_end:
        gravity.evolve_model(gravity.model_time+dt)
        if filename:
            write_set_to_file(gravity.particles, filename, 'hdf5')
        print_diagnostics(gravity.model_time, [black_hole, stars], gravity.converter) 
    gravity.stop()

def new_option_parser():
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 1,
                      help="number of diagnostics time steps [10]")
    result.add_option("-f", dest="filename", default = None,
                      help="write output filename")
    result.add_option("-t", unit=units.Myr, 
                      dest="t_end", type="float", default = 0.000001|units.Myr,
                      help="end time of the simulation [0.0000001] %unit")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

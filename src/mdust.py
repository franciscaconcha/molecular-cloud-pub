import os
import sys
import multiprocessing
import Queue
import threading
import numpy

from amuse.lab import *
from amuse.community.fractalcluster.interface import new_fractal_cluster_model
from amuse.ic.kingmodel import new_king_model
from amuse.units import units, constants

from disk_class import setup_disks_and_codes, run_disks
import FRIED_interp

G0 = 1.6e-3 * units.erg / units.s / units.cm**2


def dust_mass(star, t):
	Tm = (100. | units.K) * (star.stellar_mass.value_in(units.MSun)) ** (1. / 4.)  # Tm ~ L^1/4 ~ (M^3)^1/4
	mu = 2.33
	delta = 1e-2
	rho_g = 1. | units.g / units.cm ** 3
	a_min = 1e-8 | units.m
	t = t | units.Myr
	dt = 1000 | units.yr
	# Remove dust in a leapfrog-like integration
	# Follows the prescription of Haworth et al. 2018 (MNRAS 475)

	# Thermal speed of particles
	v_th = (8. * constants.kB * Tm / numpy.sqrt(star.disk_radius.value_in(units.AU)) / (
			numpy.pi * mu * 1.008 * constants.u)) ** (1. / 2.)
	# Disk scale height at disk edge
	Hd = (constants.kB * Tm * (1. | units.AU) ** (1. / 2.) * star.disk_radius ** (5. / 2.) / (
			mu * 1.008 * constants.u * star.stellar_mass * constants.G)) ** (1. / 2.)
	# Disk filling factor of sphere at disk edge
	F = Hd / (Hd ** 2 + star.disk_radius ** 2) ** (1. / 2.)

	star.dust_photoevap_rate = delta * star.photoevap_Mdot ** (3. / 2.) * \
	                           (v_th / (
			                           4. * numpy.pi * F * constants.G * star.stellar_mass * rho_g * a_min)) ** (
			                           1. / 2.) * \
	                           numpy.exp(-delta * (constants.G * star.stellar_mass) ** (
			                           1. / 2.) * t / (2. * star.disk_radius ** (3. / 2.)))

	# Can't entrain more dust than is available
	if star.dust_photoevap_rate > delta * star.photoevap_Mdot:
		star.dust_photoevap_rate = delta * star.photoevap_Mdot

	# Eulerian integration
	dM_dust = star.dust_photoevap_rate * dt
	if not star.disked:  # If disk is dispersed, do only half a step
		dM_dust /= 2.
	star.disk_dust_mass -= dM_dust

	# Can't have negative mass
	if star.disk_dust_mass < 0. | units.MSun:
		star.disk_dust_mass = 0. | units.MSun

	return star.disk_dust_mass


def main(open_path, N, save_path, t_end, save, nruns):
	for n in range(nruns):
		path = '{0}/{1}/disks/'.format(open_path, n)
		files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
		files = [x for x in files if '.hdf5' in x]
		files.sort(key=lambda f: float(filter(str.isdigit, f)))

		print path
		print files

		for f in files:
			stars = read_set_from_file(path + f, 'hdf5', close_file=True)
			t = float(f.split('t')[1].split('.hdf5')[0])

			#stars = stars[stars.disked]

			for s in stars[stars.disked]:
				s.disk_dust_mass = dust_mass(s, t)

			write_set_to_file(stars,
			                  '{0}/{1}/mdust/N{2}_t{3:.3f}.hdf5'.format(open_path,
			                                                      n,
			                                                      len(stars),
			                                                      t),
			                  'hdf5')


def new_option_parser():
    from amuse.units.optparse import OptionParser

    result = OptionParser()
    result.add_option("-p", dest="open_path", type="string", default='/media/fran/data1/photoevap/results',
                      help="path to results to plot [%default]")
    result.add_option("-N", dest="N", type="int", default=1000,
                      help="number of stars [%default]")
    result.add_option("-S", dest="save", type="int", default=0,
                      help="save plot? [%default]")
    result.add_option("-s", dest="save_path", type="string", default='/media/fran/data1/photoevap-paper/figures',
                      help="path to save the results [%default]")
    result.add_option("-e", dest="t_end", type="float", default='5.0',
                      help="end time to use for plots [%default]")
    result.add_option("-n", dest="nruns", type="int", default=1,
                      help="number of runs to plot for averages [%default]")
    return result

if __name__ == '__main__':
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)
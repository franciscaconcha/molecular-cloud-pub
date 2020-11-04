import os
import sys
import multiprocessing
import Queue
import threading
import numpy
#from galpy.potential import to_amuse, MWPotential2014

from amuse.lab import *
from amuse.community.fractalcluster.interface import new_fractal_cluster_model
from amuse.ic.kingmodel import new_king_model

from disk_class import setup_disks_and_codes, run_disks
import FRIED_interp

G0 = 1.6e-3 * units.erg / units.s / units.cm**2

code_queue = Queue.Queue()


def map_disk_indices_to_stars(disks):
	map = {}
	i = 0

	for d in disks:
		map[d.key] = i
		i += 1

	return map


def main(path, grid_path, ndisks, nmin, nmax):

	interpolator = FRIED_interp.FRIED_interpolator(folder=grid_path, verbosity=False)

	all_stars = read_set_from_file("{0}/gravity_stars.amuse".format(path), "hdf5", close_file=True)
	all_disked_stars = all_stars[all_stars.disked]
	disk_codes, disks = setup_disks_and_codes(all_disked_stars.key,
											  all_disked_stars.initial_disk_radius,
											  all_disked_stars.initial_disk_mass,
											  all_disked_stars.stellar_mass,
											  all_disked_stars.dispersed_mass_threshold,
											  all_disked_stars.dispersed_density_threshold,
											  ndisks,  # number of vaders
											  100,  # number of cells
											  0.05 | units.au,
											  2000 | units.au,
											  5E-3)

	all_stars.prev = False  # To keep track of born stars in timed files

	disk_indices = map_disk_indices_to_stars(disks)

	for t in times:
		born_stars = # stars born in time < t
		new_stars = # stars born in this t

		# Update disks of the newly born stars
		for s in new_stars:
			s.bright = s.stellar_mass > 1.9 | units.MSun
			s.disked = s.stellar_mass <= 1.9 | units.MSun

			if s.disked:
				s.disk_radius = 30 * (s.stellar_mass.value_in(units.MSun) ** 0.5) | units.au
				s.disk_mass = 0.1 * s.stellar_mass

				s.mass = s.stellar_mass + s.disk_mass
				s.collisional_radius = 0.02 | units.parsec
				s.encounters = 0  # Counter for dynamical encounters

				s.cumulative_truncation_mass_loss = 0.0 | units.MSun
				s.cumulative_photoevap_mass_loss = 0.0 | units.MSun

				# Saving initial G0 on small stars
				s.g0 = 0.0

				# Flag for EUV photoevaporation mass loss
				s.EUV = False

				s.dispersed_mass_threshold = 0.03 | units.MEarth  # Ansdell+2016
				s.dispersed_density_threshold = 1E-5 | units.g / units.cm ** 2  # Ingleby+ 2009
			else:
				s.disk_radius = 0.0 | units.au
				s.disk_mass = 0.0 | units.MSun

				s.mass = s.stellar_mass
				s.collisional_radius = 0.02 | units.parsec
				s.encounters = 0  # Counter for dynamical encounters

		#for s in born_stars + new_stars:
		born_stars[born_stars.disked].total_radiation = total_radiation(born_stars[born_stars.disked].key, ncores)

		# Photoevaporative mass loss in log10(MSun/yr), EUV + FUV
		born_stars[born_stars.disked].photoevap_Mdot = photoevaporation_mass_loss(born_stars[born_stars.disked].key,
																				  ncores)
		born_stars[born_stars.disked].cumulative_photoevap_mass_loss += born_stars[
																			born_stars.disked].photoevap_Mdot * dt

		# Update disks' mass loss rates before evolving them
		for k in disk_indices:
			if len(born_stars[born_stars.key == k].photoevap_Mdot) > 0:
				disks[disk_indices[k]].outer_photoevap_rate = born_stars[born_stars.key == k].photoevap_Mdot
			else:
				disks[disk_indices[k]].outer_photoevap_rate = 0.0 | units.MSun / units.yr

		# Evolve VADER disks
		# This evolution includes gas+dust evolution and external photoevaporation
		# disks_to_run = [d for d in disks if (not d.dispersed and d.born)]
		run_disks(disk_codes, [d for d in disks if not d.dispersed], dt)


		for s in born_stars + new_stars:
			# update disk parameters

		# write born_stars + new_stars to file, same as what I save in vader_cluster_parallel

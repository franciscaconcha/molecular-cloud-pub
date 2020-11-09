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


def total_radiation(indices, nc):  # indices should be list of keys of small stars
    pool = multiprocessing.Pool(processes=nc)
    total_radiation = pool.map(single_total_radiation, indices)
    pool.close()
    pool.join()
    return total_radiation


def single_total_radiation(i):
    global born_stars
    this_star = born_stars[born_stars.key == i]  # Current star to calculate total radiation on

    # Calculate the total FUV contribution of the bright stars over each small star
    total_radiation = 0.0

    if i in born_stars[born_stars.bright].key:
        print "No radiating stars yet"
        return total_radiation

    for s in born_stars[born_stars.bright]:  # For each massive/bright star
        # Calculate FUV luminosity of the bright star, in LSun
        lum = luminosity_fit(s.stellar_mass.value_in(units.MSun))

        # Calculate distance to bright star
        dist = distance(s, this_star)[0]

        # EUV regime -- Use Johnstone, Hollenbach, & Bally 1998
        dmin = 5. * 1E17 * 0.25 * numpy.sqrt(this_star.disk_radius.value_in(units.cm) / 1E14) | units.cm

        if dist < dmin:
            this_star.EUV = True
        # EUV mass loss will be calculated in mass_loss_parallel

        # FUV radiation only
        rad = lum.value_in(units.erg / units.s)
        d = dist.value_in(units.cm)
        radiation = rad / (4 * numpy.pi * d ** 2) | (units.erg / (units.s * units.cm ** 2))

        radiation_G0 = radiation.value_in(G0)
        total_radiation += radiation_G0

    return total_radiation

def distance(star1,
             star2):
    """ Return distance between star1 and star2

    :param star1: AMUSE particle
    :param star2: AMUSE particle
    :return: distance in units.parsec
    """
    return numpy.sqrt((star2.x - star1.x)**2 + (star2.y - star1.y)**2 + (star2.z - star1.z)**2)


def radiation_at_distance(rad, d):
    """ Return radiation rad at distance d

    :param rad: total radiation from star in erg/s
    :param d: distance in cm
    :return: radiation of star at distance d, in erg * s^-1 * cm^-2
    """
    return rad / (4 * numpy.pi * d**2) | (units.erg / (units.s * units.cm**2))


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


def photoevaporation_mass_loss(indices, nc):
    pool = multiprocessing.Pool(processes=nc)
    mass_losses = pool.map(single_photoevaporation_mass_loss, indices)
    pool.close()
    pool.join()
    return mass_losses


def single_photoevaporation_mass_loss(i):
    """ Return Mdot from FRIED grid.

    :param i: star key
    :return: Mdot (FUV + EUV) in MSun/yr
    """
    global born_stars, disks, disk_indices, interpolator

    this_star = born_stars[born_stars.key == i]
    this_disk = disks[disk_indices[i]]

    # FUV mass loss: interpolate from FRIED grid
    # At the start, some stars will receive no radiation because massive stars have not formed yet
    if this_star.total_radiation | G0 > 0.0 | G0:
        photoevap_Mdot_FUV = interpolator.interp_amuse(this_star.stellar_mass,
                                                       this_star.total_radiation | G0,
                                                       this_disk.disk_gas_mass,
                                                       this_disk.disk_radius)[0]

        # Check if there should be EUV mass loss as well; FRIED grid is FUV only
        # this_star.EUV is determined when the radiation on each star is calculated
        if this_star.EUV:
            # Photoevaporative mass loss in MSun/yr. Eq 20 from Johnstone, Hollenbach, & Bally 1998
            # From the paper: e ~ 3, x ~ 1.5
            photoevap_Mdot_EUV = 2. * 1E-9 * 3 * 4.12 * (this_star.disk_radius.value_in(units.cm) / 1E14) | units.MSun/units.yr
        else:
            photoevap_Mdot_EUV = 0.0 | units.MSun/units.yr

        this_star.EUV = False  # Back to false to recheck next time

        return photoevap_Mdot_FUV + photoevap_Mdot_EUV

    else:
        return 0.0 | units.MSun/units.yr


def main(open_path, grid_path, save_path, nrun, ndisks, ncores):
	"""

	:param open_path: path to results
	:param grid_path: path to fried grid
	:param ndisks: number of disks to use
	:param n: number of run to open
	:return:
	"""
	global born_stars, disks, disk_indices, interpolator

	interpolator = FRIED_interp.FRIED_interpolator(folder=grid_path, verbosity=False)

	# These are all the stars that were formed at the end of the simulation
	stars = read_set_from_file("{0}/{1}/gravity_stars.hdf5".format(open_path, nrun),
							   "hdf5", close_file=True)

	stars.born = False

	# Updating initial star parameters
	# Bright stars: no disks; emit FUV radiation
	stars[stars.stellar_mass.value_in(units.MSun) > 1.9].bright = True
	stars[stars.stellar_mass.value_in(units.MSun) > 1.9].disked = False

	# Small stars: with disks; radiation from them not considered
	stars[stars.stellar_mass.value_in(units.MSun) <= 1.9].bright = False
	stars[stars.stellar_mass.value_in(units.MSun) <= 1.9].disked = True

	stars[stars.disked].disk_radius = 30 * (stars[stars.disked].stellar_mass.value_in(units.MSun) ** 0.5) | units.au
	stars[stars.disked].disk_mass = 0.1 * stars[stars.disked].stellar_mass

	stars[stars.bright].disk_radius = 0 | units.au
	stars[stars.bright].disk_mass = 0 | units.MSun

	stars.initial_disk_radius = stars.disk_radius
	stars.initial_disk_mass = stars.disk_mass

	stars.mass = stars.stellar_mass + stars.disk_mass  # Total particle mass is stellar mass + disk mass

	# Initially all stars have the same collisional radius
	stars.collisional_radius = 0.02 | units.parsec
	stars.encounters = 0  # Counter for dynamical encounters

	stars[stars.disked].cumulative_truncation_mass_loss = 0.0 | units.MSun
	stars[stars.disked].cumulative_photoevap_mass_loss = 0.0 | units.MSun

	# Saving initial G0 on small stars
	stars.g0 = 0.0

	# Flag for EUV photoevaporation mass loss
	stars.EUV = False

	stars[stars.disked].dispersed_mass_threshold = 0.03 | units.MEarth  # Ansdell+2016
	stars[stars.disked].dispersed_density_threshold = 1E-5 | units.g / units.cm ** 2  # Ingleby+ 2009

	disk_codes, disks = setup_disks_and_codes(stars[stars.disked].key,
											  stars[stars.disked].initial_disk_radius,
											  stars[stars.disked].initial_disk_mass,
											  stars[stars.disked].stellar_mass,
											  stars[stars.disked].dispersed_mass_threshold,
											  stars[stars.disked].dispersed_density_threshold,
											  ndisks,  # number of vaders
											  100,  # number of cells
											  0.05 | units.au,
											  2000 | units.au,
											  5E-3)  # alpha

	disk_indices = map_disk_indices_to_stars(disks)

	N = len(stars)

	path = '{0}/{1}/'.format(open_path, nrun)
	files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
	files = [x for x in files if 'hydro_stars' in x]
	files.sort(key=lambda f: float(filter(str.isdigit, f)))

	t_prev = read_set_from_file("{0}/{1}".format(path, files[0]),
								"hdf5", close_file=True).get_timestamp()
	dt_prev = 0.0  # This is to take care of the last dt which can be 10 Myr bc of how I break out of the molecular cloud loop

	stellar = None
	stars[stars.bright].in_stellar = False  # This is just a workaround to not add particles twice

	for f in files:
		current_stars = read_set_from_file("{0}/{1}".format(path, f),
										   "hdf5", close_file=True)
		t = current_stars.get_timestamp()
		if t - t_prev > 5 | units.Myr:
			dt = dt_prev
		else:
			dt = t - t_prev
		print "dt = {0}".format(dt.in_(units.Myr))

		for s in current_stars:
			# I have to do this first otherwise I add particles to stellar twice
			if stars[stars.key == s.key].bright and not stars[stars.key == s.key].in_stellar:
				if stellar is None:
					stellar = SeBa()
					stellar.parameters.metallicity = 0.02
					channel_from_framework_to_stellar = stars.new_channel_to(stellar.particles)
					channel_from_stellar_to_framework = stellar.particles.new_channel_to(stars)
					stellar.particles.add_particles(stars[stars.key == s.key])
					channel_from_framework_to_stellar.copy()
					stars[stars.key == s.key].in_stellar = True
				else:
					stellar.particles.add_particles(stars[stars.key == s.key])
					channel_from_framework_to_stellar.copy()
					stars[stars.key == s.key].in_stellar = True

			stars[stars.key == s.key].born = True
			# Update positions of born stars from hydro files
			# Instead of using a dynamics code
			stars[stars.key == s.key].x = s.x
			stars[stars.key == s.key].y = s.y
			stars[stars.key == s.key].z = s.z

		# stellar evolution
		if stellar is None:
			pass
		else:
			for sp in stellar.particles:
				sp.evolve_one_step()
			channel_from_stellar_to_framework.copy()

		born_stars = stars[stars.born]
		born_stars[born_stars.disked].total_radiation = total_radiation(born_stars[born_stars.disked].key,
																	    ncores)

		# Photoevaporative mass loss in log10(MSun/yr), EUV + FUV
		born_stars[born_stars.disked].photoevap_Mdot = photoevaporation_mass_loss(born_stars[born_stars.disked].key,
																				  ncores)
		born_stars[born_stars.disked].cumulative_photoevap_mass_loss += born_stars[born_stars.disked].photoevap_Mdot * dt

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

		for this_star in born_stars[born_stars.disked]:
			# update disk parameters
			this_disk = disks[disk_indices[this_star.key]]
			this_star.disked = not this_disk.dispersed
			this_star.disk_radius = this_disk.disk_radius
			this_star.disk_dust_mass = this_disk.disk_dust_mass
			this_star.disk_gas_mass = this_disk.disk_gas_mass
			this_star.disk_mass = this_disk.disk_mass

		# write results
		write_set_to_file(born_stars,
						  '{0}/{1}/prep/N{2}_t{3:.3f}.hdf5'.format(save_path,
															       nrun,
															       N,
															       t.value_in(units.Myr)),
						  'hdf5')

		t_prev = t
		dt_prev = dt


def new_option_parser():
	from amuse.units.optparse import OptionParser
	result = OptionParser()

	# Simulation parameters
	result.add_option("-n", dest="nrun", type="int", default=0,
					  help="run number to process [%default]")
	result.add_option("-p", dest="open_path", type="string", default='.',
					  help="path to results [%default]")
	result.add_option("-s", dest="save_path", type="string", default='.',
					  help="path to save results [%default]")
	result.add_option("-f", dest="grid_path", type="string", default='data',
					  help="path for FRIED grid [%default]")
	result.add_option("-d", dest="ndisks", type="int", default=10,
					  help="number of disks [%default]")
	result.add_option("-c", dest="ncores", type="int", default=2,
					  help="number of cores [%default]")
	return result


if __name__ == '__main__':
	o, arguments = new_option_parser().parse_args()
	main(**o.__dict__)
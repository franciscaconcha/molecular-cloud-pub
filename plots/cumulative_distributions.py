import numpy
from matplotlib import pyplot
import os
import pandas

from amuse.lab import *
from amuse import io

#movie command
#"ffmpeg -framerate 5 -i {0}/%01d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {0}/movie.mp4

from legends import *
from mycolors import *


def get_cumulative_mean(data, cumulative):
	"""Returns sorted mean of the data + cumulative dist of the mean"""
	try:
		mean_data = numpy.mean(data, axis=0)
		mean_cumulative = numpy.mean(cumulative, axis=0)
		devs = numpy.std(cumulative, axis=0)
	except:
		max_len = 0
		for a in data:
			if len(a) > max_len:
				max_len = len(a)

		new_sorted = []
		for a in data:
			b = numpy.pad(a, (max_len - len(a), 0), 'constant')
			# constant_values=(min([min(r) for r in all_initial])))
			new_sorted.append(b)

		mean_data = numpy.mean(new_sorted, axis=0)

		max_len = 0
		for a in cumulative:
			if len(a) > max_len:
				max_len = len(a)

		new_sorted = []
		for a in cumulative:
			b = numpy.pad(a, (max_len - len(a), 0), 'constant')
			# constant_values=(min([min(r) for r in all_initial])))
			new_sorted.append(b)

		mean_cumulative = numpy.mean(new_sorted, axis=0)
		devs = numpy.std(new_sorted, axis=0)

	mean_data = numpy.sort(mean_data)

	cumulative_high = mean_cumulative + devs
	cumulative_low = mean_cumulative - devs

	return mean_data, mean_cumulative, cumulative_high, cumulative_low


def separate_disk_masses(open_path, save_path, t_end, N, nruns, save):
	fig1, axs1 = pyplot.subplots(1)
	fig2, axs2 = pyplot.subplots(1)

	masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	cumulative = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

	for n in range(nruns):
		path = '{0}/{1}/disks/'.format(open_path, n)
		files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
		files = [x for x in files if '.hdf5' in x]
		files.sort(key=lambda f: float(filter(str.isdigit, f)))

		first_file = files[0]
		last_file = files[-1]

		stars = read_set_from_file(path + last_file, 'hdf5', close_file=True)
		disk_masses = stars[stars.disked].disk_mass.value_in(units.MJupiter)

		sorted_masses = numpy.sort(disk_masses)
		cumulative_c = 1. * numpy.arange(len(sorted_masses)) / (len(sorted_masses) - 1)

		masses[n] = sorted_masses
		cumulative[n] = cumulative_c

	all_masses = []
	all_cumulative = []

	for n in range(nruns):
		all_masses.append(masses[n])
		all_cumulative.append(cumulative[n])

		axs1.plot(masses[n],
					cumulative[n],
					#yerr=binned_stds[n],
					color=runcolors[n],
					lw=2,
		            label="Run {0}".format(n))

		"""axs1.fill_between(local_densities[n],
							numpy.array(binned_means[n]) + numpy.array(binned_stds[n]),
							numpy.array(binned_means[n]) - numpy.array(binned_stds[n]),
							facecolor=runcolors[n],
							alpha=0.2)"""

	mean_masses, mean_cumulative, cumulative_high, cumulative_low = get_cumulative_mean(all_masses,
	                                                                                    all_cumulative
	                                                                                    )

	axs2.plot(mean_masses,
				mean_cumulative,
				color='k',
				lw=3)

	axs2.fill_between(mean_masses,
						cumulative_high,
						cumulative_low,
						facecolor='k',
						alpha=0.2)

	axs1.legend()
	axs1.set_xlabel(r'Disc mass [$\mathrm{M}_{Jup}$]')
	axs1.set_ylabel(r'$f_{\geq \mathrm{M}_{disc, dust}}$', fontsize=24)
	axs2.set_xlabel(r'Disc mass [$\mathrm{M}_{Jup}$]')
	axs2.set_ylabel(r'$f_{\geq \mathrm{M}_{disc, dust}}$', fontsize=24)

	axs1.set_xscale('log')
	axs2.set_xscale('log')


def disk_gas_masses(open_path, save_path, t_end, N, nruns, save):
	fig1, axs1 = pyplot.subplots(1)
	fig2, axs2 = pyplot.subplots(1)

	first_masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	last_masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	first_cumulative = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	last_cumulative = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

	sf_end_indices = {0: 85, 1: 86, 2: 224, 3: 127, 4: 196, 5: 171}

	for n in range(nruns):
		path = '{0}/{1}/disks/'.format(open_path, n)
		files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
		files = [x for x in files if '.hdf5' in x]
		files.sort(key=lambda f: float(filter(str.isdigit, f)))

		first_file = files[sf_end_indices[n]]
		last_file = files[-1]

		first_stars = read_set_from_file(path + first_file, 'hdf5', close_file=True)
		disk_masses = first_stars[first_stars.disked].disk_gas_mass.value_in(units.MJupiter) #/ 100
		#disk_masses = stars[stars.disked].disk_mass.value_in(units.MEarth) / 100

		sorted_masses = numpy.sort(disk_masses)
		cumulative_c = 1. * numpy.arange(len(sorted_masses)) / (len(sorted_masses) - 1)

		first_masses[n] = sorted_masses
		first_cumulative[n] = cumulative_c

		last_stars = read_set_from_file(path + last_file, 'hdf5', close_file=True)
		disk_masses = last_stars[last_stars.disked].disk_gas_mass.value_in(units.MJupiter) #/ 100
		#disk_masses = stars[stars.disked].disk_mass.value_in(units.MEarth) / 100

		sorted_masses = numpy.sort(disk_masses)
		cumulative_c = 1. * numpy.arange(len(sorted_masses)) / (len(sorted_masses) - 1)

		last_masses[n] = sorted_masses
		last_cumulative[n] = cumulative_c

	all_first_masses, all_last_masses = [], []
	all_first_cumulative, all_last_cumulative = [], []

	for n in range(nruns):
		all_first_masses.append(first_masses[n])
		all_last_masses.append(last_masses[n])
		all_first_cumulative.append(first_cumulative[n])
		all_last_cumulative.append(last_cumulative[n])

	mean_first_masses, mean_first_cumulative, high_first, low_first = get_cumulative_mean(all_first_masses, all_first_cumulative)
	mean_last_masses, mean_last_cumulative, high_last, low_last = get_cumulative_mean(all_last_masses, all_last_cumulative)

	axs2.plot(mean_first_masses,
				1.0 - mean_first_cumulative,
				c=colors['green'],
				lw=3,
				ls='--',
				label='End of star formation')

	axs2.fill_between(mean_first_masses,
						1.0 - high_first,
						1.0 - low_first,
						facecolor=colors['green'],
						alpha=0.2)

	axs2.plot(mean_last_masses,
				1.0 - mean_last_cumulative,
				c=colors['green'],
				lw=3,
				label='End of simulations')

	axs2.fill_between(mean_last_masses,
						1.0 - high_last,
						1.0 - low_last,
						facecolor=colors['green'],
						alpha=0.2)

	axs1.legend()
	axs2.legend()

	#axs2.set_xlim([0.04, 500])
	axs2.set_ylim([0.0, 1.0])

	axs1.set_xlabel(r'Disc gas mass [$\mathrm{M}_{Jup}$]')
	axs1.set_ylabel(r'$f_{\geq \mathrm{M}_{disc, gas}}$', fontsize=24)
	axs2.set_xlabel(r'Disc gas mass [$\mathrm{M}_{Jup}$]')
	axs2.set_ylabel(r'$f_{\geq \mathrm{M}_{disc, gas}}$', fontsize=24)

	axs1.set_xscale('log')
	axs2.set_xscale('log')

	# Plotting observational data
    # OMC-2
    # Masses are sorted in the table
	"""OMC2 = pandas.read_csv('data/OMC-2_vanTerwisga2019.txt',
                           sep='\xc2\xb1',
                           names=['disk_mass', 'error'],
                           skiprows=3,
                           dtype=numpy.float64)

	OMC2_higher, OMC2_lower = [], []

	for m, e in zip(OMC2.disk_mass, OMC2.error):
		OMC2_higher.append(m + e)
		OMC2_lower.append(m - e)

	OMC2_cumulative = 1. * numpy.arange(len(OMC2.disk_mass)) / (len(OMC2.disk_mass) - 1)
	OMC2_error_cumulative = 1. * numpy.arange(len(OMC2.error)) / (len(OMC2.error) - 1)

	axs2.plot(OMC2.disk_mass, OMC2_cumulative, c='r', lw=3, label="OMC-2")

	axs2.fill_betweenx(OMC2_error_cumulative,
                      OMC2_higher,
                      OMC2_lower,
                      facecolor='r',
                      alpha=0.2)

	# ONC
	ONC_Eisner2018 = pandas.read_csv('data/ONC_Eisner2018.txt',
									sep='&',
									names=['ID', 'alpha', 'delta', 'M_star', 'F_{\rm \lambda 850 \mu m}', 'F_{\rm dust}',
	                                        'M_dust', 'R_disk'],
									skiprows=4)

	ONC_masses, ONC_error = [], []

	for me in ONC_Eisner2018.M_dust:
		m, e = me.split('$\pm$')
		ONC_masses.append(float(m))
		ONC_error.append(float(e))

	ONC_Mann2014 = pandas.read_csv('data/ONC_Mann2014.txt',
	                                 sep='\t',
	                                 names=['Field', 'Name', 'alpha', 'delta', 'M_star',
	                                        'F_{\rm \lambda 850 \mu m}', 'F_{\rm dust}',
	                                        'M_dust', 'd', 'Maj', 'Min', 'P.A.', 'Notes'],
	                                 skiprows=7)

	for me in ONC_Mann2014.M_dust:
		try:
			m, e = me.split('+or-')
		except ValueError:  # For the *one* row that is different
			m = me.split('<or=')[1]
			e = 0.0

		ONC_masses.append(float(m))
		ONC_error.append(float(e))

	ONC_higher = numpy.array(ONC_masses) + numpy.array(ONC_error)
	ONC_lower = numpy.array(ONC_masses) - numpy.array(ONC_error)

	ONC_masses.sort()

	ONC_masses_cumulative = 1. * numpy.arange(len(ONC_masses)) / (len(ONC_masses) - 1)
	ONC_error_cumulative = 1. * numpy.arange(len(ONC_error)) / (len(ONC_error) - 1)

	axs2.plot(ONC_masses[::-1], ONC_masses_cumulative, c='b', lw=3, label="ONC")
	axs2.plot(ONC_masses[::-1], ONC_higher, c='k', lw=3, label="ONC")

	axs2.fill_between(ONC_error_cumulative,
	                  ONC_higher,
	                  ONC_lower,
	                  facecolor='b',
	                  alpha=0.2)

	# Lupus
	Lupus_Ansdell2016 = pandas.read_csv('data/Lupus_Ansdell2016.txt',
	                                 sep='&',
	                                 names=['Name', 'RAh', 'DE-',
	                                        'Fcont', 'e_Fcont',
	                                        'rms', 'a', 'e_a', 'PosAng', 'e_PosAng', 'i', 'e_i',
	                                        'M_dust', 'e_M_dust'],
	                                 skiprows=31)

	Lupus_Ansdell2018 = pandas.read_csv('data/Lupus_Ansdell2018.txt',
	                                 sep='&',
	                                 names=['everything_else',
	                                        'M_dust', 'e_M_dust'],
	                                 skiprows=47)

	Lupus_masses = numpy.concatenate([Lupus_Ansdell2016['M_dust'].to_numpy(), Lupus_Ansdell2018['M_dust'].to_numpy()])
	Lupus_errors = numpy.concatenate([Lupus_Ansdell2016['e_M_dust'].to_numpy(), Lupus_Ansdell2018['e_M_dust'].to_numpy()])

	Lupus_sorted_errors = [x for _, x in sorted(zip(Lupus_masses, Lupus_errors))]
	Lupus_masses.sort()

	Lupus_higher, Lupus_lower = [], []

	for m, e in zip(Lupus_masses, Lupus_sorted_errors):
		Lupus_higher.append(m + e)
		Lupus_lower.append(m - e)

	#Lupus_higher = Lupus_masses + Lupus_sorted_errors
	#Lupus_lower = Lupus_masses - Lupus_sorted_errors

	Lupus_masses_cumulative = 1. * numpy.arange(len(Lupus_masses)) / (len(Lupus_masses) - 1)
	Lupus_errors_cumulative = 1. * numpy.arange(len(Lupus_sorted_errors)) / (len(Lupus_sorted_errors) - 1)

	axs2.plot(Lupus_masses[::-1], Lupus_masses_cumulative, c='g', lw=3, label="Lupus")
	axs2.plot(Lupus_masses[::-1], Lupus_higher[::-1], c='k', lw=3, label="Lupus")

	axs2.fill_between(Lupus_errors_cumulative,
	                  Lupus_higher,
	                  Lupus_lower,
	                  facecolor='k',
	                  alpha=0.2)"""

	if save:
		pyplot.savefig('{0}/cumulative_dust_mass.png'.format(save_path))


def disk_dust_masses(open_path, save_path, t_end, N, nruns, save):
	fig1, axs1 = pyplot.subplots(1)
	fig2, axs2 = pyplot.subplots(1)

	first_masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	last_masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	first_cumulative = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	last_cumulative = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

	sf_end_indices = {0: 85, 1: 86, 2: 224, 3: 127, 4: 196, 5: 171}

	for n in range(nruns):
		path = '{0}/{1}/disks/'.format(open_path, n)
		files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
		files = [x for x in files if '.hdf5' in x]
		files.sort(key=lambda f: float(filter(str.isdigit, f)))

		first_file = files[sf_end_indices[n]]
		last_file = files[-1]

		first_stars = read_set_from_file(path + first_file, 'hdf5', close_file=True)
		disk_masses = first_stars[first_stars.disked].disk_dust_mass.value_in(units.MEarth) #/ 100
		#disk_masses = stars[stars.disked].disk_mass.value_in(units.MEarth) / 100

		sorted_masses = numpy.sort(disk_masses)
		cumulative_c = 1. * numpy.arange(len(sorted_masses)) / (len(sorted_masses) - 1)

		first_masses[n] = sorted_masses
		first_cumulative[n] = cumulative_c

		last_stars = read_set_from_file(path + last_file, 'hdf5', close_file=True)
		disk_masses = last_stars[last_stars.disked].disk_dust_mass.value_in(units.MEarth) #/ 100
		#disk_masses = stars[stars.disked].disk_mass.value_in(units.MEarth) / 100

		sorted_masses = numpy.sort(disk_masses)
		cumulative_c = 1. * numpy.arange(len(sorted_masses)) / (len(sorted_masses) - 1)

		last_masses[n] = sorted_masses
		last_cumulative[n] = cumulative_c

	all_first_masses, all_last_masses = [], []
	all_first_cumulative, all_last_cumulative = [], []

	for n in range(nruns):
		all_first_masses.append(first_masses[n])
		all_last_masses.append(last_masses[n])
		all_first_cumulative.append(first_cumulative[n])
		all_last_cumulative.append(last_cumulative[n])

	mean_first_masses, mean_first_cumulative, high_first, low_first = get_cumulative_mean(all_first_masses, all_first_cumulative)
	mean_last_masses, mean_last_cumulative, high_last, low_last = get_cumulative_mean(all_last_masses, all_last_cumulative)

	axs2.plot(mean_first_masses,
				1.0 - mean_first_cumulative,
				c=colors['navy'],
				lw=3,
				ls='--',
				label='End of star formation')

	axs2.fill_between(mean_first_masses,
						1.0 - high_first,
						1.0 - low_first,
						facecolor=colors['navy'],
						alpha=0.2)

	axs2.plot(mean_last_masses,
				1.0 - mean_last_cumulative,
				c=colors['navy'],
				lw=3,
				label='End of simulations')

	axs2.fill_between(mean_last_masses,
						1.0 - high_last,
						1.0 - low_last,
						facecolor=colors['navy'],
						alpha=0.2)

	axs1.legend()
	axs2.legend()

	#axs2.set_xlim([0.04, 500])
	axs2.set_ylim([0.0, 1.0])

	axs1.set_xlabel(r'Disc dust mass [$\mathrm{M}_{\oplus}$]')
	axs1.set_ylabel(r'$f_{\geq \mathrm{M}_{disc, dust}}$', fontsize=24)
	axs2.set_xlabel(r'Disc dust mass [$\mathrm{M}_{\oplus}$]')
	axs2.set_ylabel(r'$f_{\geq \mathrm{M}_{disc, dust}}$', fontsize=24)

	axs1.set_xscale('log')
	axs2.set_xscale('log')

	# Plotting observational data
    # OMC-2
    # Masses are sorted in the table
	"""OMC2 = pandas.read_csv('data/OMC-2_vanTerwisga2019.txt',
                           sep='\xc2\xb1',
                           names=['disk_mass', 'error'],
                           skiprows=3,
                           dtype=numpy.float64)

	OMC2_higher, OMC2_lower = [], []

	for m, e in zip(OMC2.disk_mass, OMC2.error):
		OMC2_higher.append(m + e)
		OMC2_lower.append(m - e)

	OMC2_cumulative = 1. * numpy.arange(len(OMC2.disk_mass)) / (len(OMC2.disk_mass) - 1)
	OMC2_error_cumulative = 1. * numpy.arange(len(OMC2.error)) / (len(OMC2.error) - 1)

	axs2.plot(OMC2.disk_mass, OMC2_cumulative, c='r', lw=3, label="OMC-2")

	axs2.fill_betweenx(OMC2_error_cumulative,
                      OMC2_higher,
                      OMC2_lower,
                      facecolor='r',
                      alpha=0.2)

	# ONC
	ONC_Eisner2018 = pandas.read_csv('data/ONC_Eisner2018.txt',
									sep='&',
									names=['ID', 'alpha', 'delta', 'M_star', 'F_{\rm \lambda 850 \mu m}', 'F_{\rm dust}',
	                                        'M_dust', 'R_disk'],
									skiprows=4)

	ONC_masses, ONC_error = [], []

	for me in ONC_Eisner2018.M_dust:
		m, e = me.split('$\pm$')
		ONC_masses.append(float(m))
		ONC_error.append(float(e))

	ONC_Mann2014 = pandas.read_csv('data/ONC_Mann2014.txt',
	                                 sep='\t',
	                                 names=['Field', 'Name', 'alpha', 'delta', 'M_star',
	                                        'F_{\rm \lambda 850 \mu m}', 'F_{\rm dust}',
	                                        'M_dust', 'd', 'Maj', 'Min', 'P.A.', 'Notes'],
	                                 skiprows=7)

	for me in ONC_Mann2014.M_dust:
		try:
			m, e = me.split('+or-')
		except ValueError:  # For the *one* row that is different
			m = me.split('<or=')[1]
			e = 0.0

		ONC_masses.append(float(m))
		ONC_error.append(float(e))

	ONC_higher = numpy.array(ONC_masses) + numpy.array(ONC_error)
	ONC_lower = numpy.array(ONC_masses) - numpy.array(ONC_error)

	ONC_masses.sort()

	ONC_masses_cumulative = 1. * numpy.arange(len(ONC_masses)) / (len(ONC_masses) - 1)
	ONC_error_cumulative = 1. * numpy.arange(len(ONC_error)) / (len(ONC_error) - 1)

	axs2.plot(ONC_masses[::-1], ONC_masses_cumulative, c='b', lw=3, label="ONC")
	axs2.plot(ONC_masses[::-1], ONC_higher, c='k', lw=3, label="ONC")

	axs2.fill_between(ONC_error_cumulative,
	                  ONC_higher,
	                  ONC_lower,
	                  facecolor='b',
	                  alpha=0.2)

	# Lupus
	Lupus_Ansdell2016 = pandas.read_csv('data/Lupus_Ansdell2016.txt',
	                                 sep='&',
	                                 names=['Name', 'RAh', 'DE-',
	                                        'Fcont', 'e_Fcont',
	                                        'rms', 'a', 'e_a', 'PosAng', 'e_PosAng', 'i', 'e_i',
	                                        'M_dust', 'e_M_dust'],
	                                 skiprows=31)

	Lupus_Ansdell2018 = pandas.read_csv('data/Lupus_Ansdell2018.txt',
	                                 sep='&',
	                                 names=['everything_else',
	                                        'M_dust', 'e_M_dust'],
	                                 skiprows=47)

	Lupus_masses = numpy.concatenate([Lupus_Ansdell2016['M_dust'].to_numpy(), Lupus_Ansdell2018['M_dust'].to_numpy()])
	Lupus_errors = numpy.concatenate([Lupus_Ansdell2016['e_M_dust'].to_numpy(), Lupus_Ansdell2018['e_M_dust'].to_numpy()])

	Lupus_sorted_errors = [x for _, x in sorted(zip(Lupus_masses, Lupus_errors))]
	Lupus_masses.sort()

	Lupus_higher, Lupus_lower = [], []

	for m, e in zip(Lupus_masses, Lupus_sorted_errors):
		Lupus_higher.append(m + e)
		Lupus_lower.append(m - e)

	#Lupus_higher = Lupus_masses + Lupus_sorted_errors
	#Lupus_lower = Lupus_masses - Lupus_sorted_errors

	Lupus_masses_cumulative = 1. * numpy.arange(len(Lupus_masses)) / (len(Lupus_masses) - 1)
	Lupus_errors_cumulative = 1. * numpy.arange(len(Lupus_sorted_errors)) / (len(Lupus_sorted_errors) - 1)

	axs2.plot(Lupus_masses[::-1], Lupus_masses_cumulative, c='g', lw=3, label="Lupus")
	axs2.plot(Lupus_masses[::-1], Lupus_higher[::-1], c='k', lw=3, label="Lupus")

	axs2.fill_between(Lupus_errors_cumulative,
	                  Lupus_higher,
	                  Lupus_lower,
	                  facecolor='k',
	                  alpha=0.2)"""

	if save:
		pyplot.savefig('{0}/cumulative_dust_mass.png'.format(save_path))


def disk_masses(open_path, save_path, t_end, N, nruns, save):
	fig2, axs2 = pyplot.subplots(1)

	first_masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	last_masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	first_cumulative = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	last_cumulative = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

	sf_end_indices = {0: 85, 1: 86, 2: 224, 3: 127, 4: 196, 5: 171}

	for n in range(nruns):
		path = '{0}/{1}/disks/'.format(open_path, n)
		files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
		files = [x for x in files if '.hdf5' in x]
		files.sort(key=lambda f: float(filter(str.isdigit, f)))

		first_file = files[sf_end_indices[n]]
		last_file = files[-1]

		first_stars = read_set_from_file(path + first_file, 'hdf5', close_file=True)
		disk_masses = first_stars[first_stars.disked].disk_dust_mass.value_in(units.MEarth) #/ 100
		#disk_masses = stars[stars.disked].disk_mass.value_in(units.MEarth) / 100

		sorted_masses = numpy.sort(disk_masses)
		cumulative_c = 1. * numpy.arange(len(sorted_masses)) / (len(sorted_masses) - 1)

		first_masses[n] = sorted_masses
		first_cumulative[n] = cumulative_c

		last_stars = read_set_from_file(path + last_file, 'hdf5', close_file=True)
		disk_masses = last_stars[last_stars.disked].disk_dust_mass.value_in(units.MEarth) #/ 100
		#disk_masses = stars[stars.disked].disk_mass.value_in(units.MEarth) / 100

		sorted_masses = numpy.sort(disk_masses)
		cumulative_c = 1. * numpy.arange(len(sorted_masses)) / (len(sorted_masses) - 1)

		last_masses[n] = sorted_masses
		last_cumulative[n] = cumulative_c

	all_first_masses, all_last_masses = [], []
	all_first_cumulative, all_last_cumulative = [], []

	for n in range(nruns):
		all_first_masses.append(first_masses[n])
		all_last_masses.append(last_masses[n])
		all_first_cumulative.append(first_cumulative[n])
		all_last_cumulative.append(last_cumulative[n])

		"""axs1.plot(masses[n],
					cumulative[n],
					#yerr=binned_stds[n],
					color=runcolors[n],
					lw=3,
		            label="Run {0}".format(n))"""

		"""axs1.fill_between(local_densities[n],
							numpy.array(binned_means[n]) + numpy.array(binned_stds[n]),
							numpy.array(binned_means[n]) - numpy.array(binned_stds[n]),
							facecolor=runcolors[n],
							alpha=0.2)"""

	mean_first_masses, mean_first_cumulative, high_first, low_first = get_cumulative_mean(all_first_masses, all_first_cumulative)
	mean_last_masses, mean_last_cumulative, high_last, low_last = get_cumulative_mean(all_last_masses, all_last_cumulative)

	axs2.plot(mean_first_masses,
				1.0 - mean_first_cumulative,
				c=colors['navy'],
				lw=3,
				ls='--',
				label='End of star formation')

	axs2.fill_between(mean_first_masses,
						1.0 - high_first,
						1.0 - low_first,
						facecolor=colors['navy'],
						alpha=0.2)

	axs2.plot(mean_last_masses,
				1.0 - mean_last_cumulative,
				c=colors['navy'],
				lw=3,
				label='End of simulations')

	axs2.fill_between(mean_last_masses,
						1.0 - high_last,
						1.0 - low_last,
						facecolor=colors['navy'],
						alpha=0.2)

	first_masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	last_masses = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	first_cumulative = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	last_cumulative = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

	sf_end_indices = {0: 85, 1: 86, 2: 224, 3: 127, 4: 196, 5: 171}

	for n in range(nruns):
		path = '{0}/{1}/disks/'.format(open_path, n)
		files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
		files = [x for x in files if '.hdf5' in x]
		files.sort(key=lambda f: float(filter(str.isdigit, f)))

		first_file = files[sf_end_indices[n]]
		last_file = files[-1]

		first_stars = read_set_from_file(path + first_file, 'hdf5', close_file=True)
		disk_masses = first_stars[first_stars.disked].disk_gas_mass.value_in(units.MEarth)  # / 100
		# disk_masses = stars[stars.disked].disk_mass.value_in(units.MEarth) / 100

		sorted_masses = numpy.sort(disk_masses)
		cumulative_c = 1. * numpy.arange(len(sorted_masses)) / (len(sorted_masses) - 1)

		first_masses[n] = sorted_masses
		first_cumulative[n] = cumulative_c

		last_stars = read_set_from_file(path + last_file, 'hdf5', close_file=True)
		disk_masses = last_stars[last_stars.disked].disk_gas_mass.value_in(units.MEarth)  # / 100
		# disk_masses = stars[stars.disked].disk_mass.value_in(units.MEarth) / 100

		sorted_masses = numpy.sort(disk_masses)
		cumulative_c = 1. * numpy.arange(len(sorted_masses)) / (len(sorted_masses) - 1)

		last_masses[n] = sorted_masses
		last_cumulative[n] = cumulative_c

	all_first_masses, all_last_masses = [], []
	all_first_cumulative, all_last_cumulative = [], []

	for n in range(nruns):
		all_first_masses.append(first_masses[n])
		all_last_masses.append(last_masses[n])
		all_first_cumulative.append(first_cumulative[n])
		all_last_cumulative.append(last_cumulative[n])

	mean_first_masses, mean_first_cumulative, high_first, low_first = get_cumulative_mean(all_first_masses,
																							all_first_cumulative)
	mean_last_masses, mean_last_cumulative, high_last, low_last = get_cumulative_mean(all_last_masses,
																						all_last_cumulative)

	axs2.plot(mean_first_masses,
				1.0 - mean_first_cumulative,
				c=colors['turquoise'],
				lw=3,
				ls='--',
				label='End of star formation')

	axs2.fill_between(mean_first_masses,
						1.0 - high_first,
						1.0 - low_first,
						facecolor=colors['turquoise'],
						alpha=0.2)

	axs2.plot(mean_last_masses,
				1.0 - mean_last_cumulative,
				c=colors['turquoise'],
				lw=3,
				label='End of simulations')

	axs2.fill_between(mean_last_masses,
						1.0 - high_last,
						1.0 - low_last,
						facecolor=colors['turquoise'],
						alpha=0.2)

	first_legend = axs2.legend([GasSolidDashedObject(), DustSolidDashedObject()],
								['Gas mass', 'Dust mass'],
								handler_map={GasSolidDashedObject: GasSolidDashedObjectHandler(),
											DustSolidDashedObject: DustSolidDashedObjectHandler(),},
								loc='lower left',
								bbox_to_anchor=(0.004, 0.18),
								#ncol=2,
								fontsize=22, framealpha=0.4)

	axs2.add_artist(first_legend)

	axs2.legend([DashedObject(), SolidObject()],
				['End of star formation', 'End of simulations'],
				handler_map={DashedObject: DashedObjectHandler(),
							SolidObject: SolidObjectHandler(), },
				loc='lower left',
				#bbox_to_anchor=(0.52, 0.2),
				# ncol=2,
				fontsize=22, framealpha=0.4)

	# axs2.set_xlim([0.04, 500])
	axs2.set_ylim([0.0, 1.0])

	axs2.set_xlim([0.01, 50000])
	axs2.set_ylim([0.0, 1.0])

	axs2.set_xlabel(r'Disc mass [$\mathrm{M}_{\oplus}$]')
	axs2.set_ylabel(r'$f_{\geq \mathrm{M}_{disc}}$', fontsize=24)

	axs2.set_xscale('log')

	# Plotting observational data
    # OMC-2
    # Masses are sorted in the table
	"""OMC2 = pandas.read_csv('data/OMC-2_vanTerwisga2019.txt',
                           sep='\xc2\xb1',
                           names=['disk_mass', 'error'],
                           skiprows=3,
                           dtype=numpy.float64)

	OMC2_higher, OMC2_lower = [], []

	for m, e in zip(OMC2.disk_mass, OMC2.error):
		OMC2_higher.append(m + e)
		OMC2_lower.append(m - e)

	OMC2_cumulative = 1. * numpy.arange(len(OMC2.disk_mass)) / (len(OMC2.disk_mass) - 1)
	OMC2_error_cumulative = 1. * numpy.arange(len(OMC2.error)) / (len(OMC2.error) - 1)

	axs2.plot(OMC2.disk_mass, OMC2_cumulative, c='r', lw=3, label="OMC-2")

	axs2.fill_betweenx(OMC2_error_cumulative,
                      OMC2_higher,
                      OMC2_lower,
                      facecolor='r',
                      alpha=0.2)

	# ONC
	ONC_Eisner2018 = pandas.read_csv('data/ONC_Eisner2018.txt',
									sep='&',
									names=['ID', 'alpha', 'delta', 'M_star', 'F_{\rm \lambda 850 \mu m}', 'F_{\rm dust}',
	                                        'M_dust', 'R_disk'],
									skiprows=4)

	ONC_masses, ONC_error = [], []

	for me in ONC_Eisner2018.M_dust:
		m, e = me.split('$\pm$')
		ONC_masses.append(float(m))
		ONC_error.append(float(e))

	ONC_Mann2014 = pandas.read_csv('data/ONC_Mann2014.txt',
	                                 sep='\t',
	                                 names=['Field', 'Name', 'alpha', 'delta', 'M_star',
	                                        'F_{\rm \lambda 850 \mu m}', 'F_{\rm dust}',
	                                        'M_dust', 'd', 'Maj', 'Min', 'P.A.', 'Notes'],
	                                 skiprows=7)

	for me in ONC_Mann2014.M_dust:
		try:
			m, e = me.split('+or-')
		except ValueError:  # For the *one* row that is different
			m = me.split('<or=')[1]
			e = 0.0

		ONC_masses.append(float(m))
		ONC_error.append(float(e))

	ONC_higher = numpy.array(ONC_masses) + numpy.array(ONC_error)
	ONC_lower = numpy.array(ONC_masses) - numpy.array(ONC_error)

	ONC_masses.sort()

	ONC_masses_cumulative = 1. * numpy.arange(len(ONC_masses)) / (len(ONC_masses) - 1)
	ONC_error_cumulative = 1. * numpy.arange(len(ONC_error)) / (len(ONC_error) - 1)

	axs2.plot(ONC_masses[::-1], ONC_masses_cumulative, c='b', lw=3, label="ONC")
	axs2.plot(ONC_masses[::-1], ONC_higher, c='k', lw=3, label="ONC")

	axs2.fill_between(ONC_error_cumulative,
	                  ONC_higher,
	                  ONC_lower,
	                  facecolor='b',
	                  alpha=0.2)

	# Lupus
	Lupus_Ansdell2016 = pandas.read_csv('data/Lupus_Ansdell2016.txt',
	                                 sep='&',
	                                 names=['Name', 'RAh', 'DE-',
	                                        'Fcont', 'e_Fcont',
	                                        'rms', 'a', 'e_a', 'PosAng', 'e_PosAng', 'i', 'e_i',
	                                        'M_dust', 'e_M_dust'],
	                                 skiprows=31)

	Lupus_Ansdell2018 = pandas.read_csv('data/Lupus_Ansdell2018.txt',
	                                 sep='&',
	                                 names=['everything_else',
	                                        'M_dust', 'e_M_dust'],
	                                 skiprows=47)

	Lupus_masses = numpy.concatenate([Lupus_Ansdell2016['M_dust'].to_numpy(), Lupus_Ansdell2018['M_dust'].to_numpy()])
	Lupus_errors = numpy.concatenate([Lupus_Ansdell2016['e_M_dust'].to_numpy(), Lupus_Ansdell2018['e_M_dust'].to_numpy()])

	Lupus_sorted_errors = [x for _, x in sorted(zip(Lupus_masses, Lupus_errors))]
	Lupus_masses.sort()

	Lupus_higher, Lupus_lower = [], []

	for m, e in zip(Lupus_masses, Lupus_sorted_errors):
		Lupus_higher.append(m + e)
		Lupus_lower.append(m - e)

	#Lupus_higher = Lupus_masses + Lupus_sorted_errors
	#Lupus_lower = Lupus_masses - Lupus_sorted_errors

	Lupus_masses_cumulative = 1. * numpy.arange(len(Lupus_masses)) / (len(Lupus_masses) - 1)
	Lupus_errors_cumulative = 1. * numpy.arange(len(Lupus_sorted_errors)) / (len(Lupus_sorted_errors) - 1)

	axs2.plot(Lupus_masses[::-1], Lupus_masses_cumulative, c='g', lw=3, label="Lupus")
	axs2.plot(Lupus_masses[::-1], Lupus_higher[::-1], c='k', lw=3, label="Lupus")

	axs2.fill_between(Lupus_errors_cumulative,
	                  Lupus_higher,
	                  Lupus_lower,
	                  facecolor='k',
	                  alpha=0.2)"""

	if save:
		pyplot.savefig('{0}/cumulative_disk_mass.png'.format(save_path))


def main(open_path, N, save_path, t_end, save, nruns):
    # My own stylesheet, comment out if not needed
    pyplot.style.use('paper')

    #disk_dust_masses(open_path, save_path, t_end, N, nruns, save)
    #disk_gas_masses(open_path, save_path, t_end, N, nruns, save)
    disk_masses(open_path, save_path, t_end, N, nruns, save)

    if not save:
        pyplot.show()


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
    result.add_option("-n", dest="nruns", type="int", default=5,
                      help="number of runs to plot for averages [%default]")
    return result


if __name__ == '__main__':
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)

import numpy
from matplotlib import pyplot
from scipy import stats
from sklearn.neighbors import KDTree
import os

from amuse.lab import *
from amuse import io

#movie command
#"ffmpeg -framerate 5 -i {0}/%01d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {0}/movie.mp4

from legends import *  # My own custom legend definitions
from mycolors import *


def mass_vs_local_density(open_path, save_path, t_end, N, nruns, save):
	fig1, axs1 = pyplot.subplots(1)
	fig2, axs2 = pyplot.subplots(1)

	times = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	binned_means = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	binned_stds = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
	local_densities = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}

	for n in range(nruns):
		path = '{0}/{1}/disks/'.format(open_path, n)
		files = os.listdir(path)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
		files = [x for x in files if '.hdf5' in x]
		files.sort(key=lambda f: float(filter(str.isdigit, f)))

		first_file = files[0]
		last_file = files[-1]

		stars = read_set_from_file(path + last_file, 'hdf5', close_file=True)
		t = float(last_file.split('t')[1].split('.hdf5')[0])
		times[n].append(t)

		stars = stars[stars.disked]

		positions = []
		for s in stars:
			positions.append(numpy.array([s.x.value_in(units.parsec),
										  s.y.value_in(units.parsec),#]))#,
										  s.z.value_in(units.parsec)]))

		positions = numpy.array(positions)

		tree = KDTree(positions)
		nearest_dist, nearest_ind = tree.query(positions, k=5)
		# print(nearest_dist)  # drop id; assumes sorted -> see args!
		# print(nearest_ind)

		for i in range(len(stars)):
			s = stars[i]
			distances_to_neighbours = nearest_dist[i]
			max_distance = max(distances_to_neighbours)
			s.local_density = 5. / ((4. / 3.) * numpy.pi * max_distance ** 3)

		#disk_masses = stars.disk_dust_mass.value_in(units.MEarth)# / 100
		disk_masses = stars.disk_mass.value_in(units.MEarth)# / 100
		masses_sorted_by_local_dens = [float(x) for _, x in sorted(zip(stars.local_density, disk_masses))]
		sorted_local_dens = sorted(stars.local_density)

		d = 100  # number of points in each bin
		for i in range(len(masses_sorted_by_local_dens)):
			if len(masses_sorted_by_local_dens) - (i + d) > 0:
				binned_means[n].append(numpy.mean(masses_sorted_by_local_dens[i:i+d]))
				binned_stds[n].append(stats.sem(masses_sorted_by_local_dens[i:i+d]))
				local_densities[n].append(numpy.mean(sorted_local_dens[i:i+d]))
			else:
				#print "end"
				binned_means[n].append(numpy.mean(masses_sorted_by_local_dens[i:]))
				binned_stds[n].append(stats.sem(masses_sorted_by_local_dens[i:]))
				local_densities[n].append(numpy.mean(sorted_local_dens[i:i+d]))
				break

	all_binned_means = []
	all_stds = []
	all_local_densities = []

	for n in range(nruns):
		all_binned_means.append(binned_means[n])
		all_stds.append(binned_stds[n])
		all_local_densities.append(local_densities[n])

		axs1.plot(local_densities[n],
					binned_means[n],
					#yerr=binned_stds[n],
					color=runcolors[n],
					lw=2)

		axs1.fill_between(local_densities[n],
							numpy.array(binned_means[n]) + numpy.array(binned_stds[n]),
							numpy.array(binned_means[n]) - numpy.array(binned_stds[n]),
							facecolor=runcolors[n],
							alpha=0.2)

	try:
		all_means = numpy.mean(all_binned_means, axis=0)
		devs = numpy.mean(all_stds, axis=0)
	except:
		max_len = 0
		for a in all_binned_means:
			if len(a) > max_len:
				max_len = len(a)

		new_sorted = []
		new_stds = []
		for a in all_binned_means:
			b = numpy.pad(a, (max_len - len(a), 0), 'constant')
			# constant_values=(min([min(r) for r in all_initial])))
			new_sorted.append(b)
		for a in all_stds:
			b = numpy.pad(a, (max_len - len(a), 0), 'constant')
			# constant_values=(min([min(r) for r in all_initial])))
			new_stds.append(b)

		all_means = numpy.mean(new_sorted, axis=0)
		devs = numpy.mean(new_stds, axis=0)

	all_means_high = all_means + devs
	all_means_low = all_means - devs

	try:
		locdens_means = numpy.mean(all_local_densities, axis=0)
	except:
		max_len = 0
		for a in all_local_densities:
			if len(a) > max_len:
				max_len = len(a)

		new_sorted = []
		for a in all_local_densities:
			b = numpy.pad(a, (max_len - len(a), 0), 'constant')
			# constant_values=(min([min(r) for r in all_initial])))
			new_sorted.append(b)
		locdens_means = numpy.mean(new_sorted, axis=0)

	axs2.plot(locdens_means,
				all_means,
				color='k',
				lw=3)

	axs2.fill_between(locdens_means,
						all_means_high,
						all_means_low,
						facecolor='k',
						alpha=0.2)


	axs1.set_xlabel(r'Local stellar number density [pc$^{-3}$]')
	axs1.set_ylabel(r'Binned mean disc dust mass [$\mathrm{M}_{\oplus}$]')

	axs1.set_xscale('log')
	axs1.set_yscale('log')

	axs2.set_xlabel(r'Local stellar number density [pc$^{-3}$]')
	axs2.set_ylabel(r'Binned mean disc dust mass [$\mathrm{M}_{\oplus}$]')

	axs2.set_xscale('log')
	axs2.set_yscale('log')

	if save:
		pyplot.savefig('{0}/2D_dustmass_localdensity.png'.format(save_path))


def main(open_path, N, save_path, t_end, save, nruns, all_in_one):

	# My own stylesheet, comment out if not needed
	pyplot.style.use('paper')

	mass_vs_local_density(open_path, save_path, t_end, N, nruns, save)

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
	result.add_option("-n", dest="nruns", type="int", default=1,
					  help="number of runs to plot for averages [%default]")
	result.add_option("-a", dest="all_in_one", type="int", default=1,
					  help="plot all in one single figure [%default]")
	return result


if __name__ == '__main__':
	o, arguments = new_option_parser().parse_args()
	main(**o.__dict__)

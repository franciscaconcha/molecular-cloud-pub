import numpy
from matplotlib import pyplot
import os

from amuse.lab import *


def mean_sink_size_vs_Nsph(path, save_path, Mcloud, Rcloud):

    Mcloud = [int(Mcloud.value_in(units.MSun))]#, 7500, 15000]
    SFE = [40, 25, 10]

    for M in Mcloud:

        Nsph = [4000, 8000, 16000, 32000]
        Nruns = [1, 2, 3, 4, 5, 6, 7, 8]

        plots = []

        for N in Nsph:
            sizes = []
            for r in Nruns:
                if N * r <= 32000:
                    filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                                      M,
                                                                      int(Rcloud.value_in(units.parsec)),
                                                                      N,
                                                                      r)

                    files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'

                    for f in files:
                        if 'sink' in f:
                            sink_particles = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)
                            sizes.append(numpy.mean(sink_particles.radius.value_in(units.parsec)))
            plots.append(numpy.mean(sizes))

        pyplot.plot(Nsph, plots, lw=2, label='SFE={0}\%'.format(SFE[Mcloud.index(M)]))

    pyplot.xlabel(r'$N_\mathrm{SPH}$')
    pyplot.ylabel(r'$<R_\mathrm{sink}>$ [pc]')
    pyplot.xticks(Nsph)
    pyplot.legend()
    pyplot.savefig('{0}/sink_size.png'.format(save_path))
    pyplot.show()


def mean_sink_mass_vs_Nsph(path, save_path, Mcloud, Rcloud):

    Mcloud = [int(Mcloud.value_in(units.MSun))]#, 7500, 15000]
    SFE = [40, 25, 10]

    for M in Mcloud:

        Nsph = [4000, 8000, 16000, 32000]
        Nruns = [1, 2, 3, 4, 5, 6, 7, 8]

        plots = []

        for N in Nsph:
            masses = []
            for r in Nruns:
                if N * r <= 32000:
                    filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                                      M,
                                                                      int(Rcloud.value_in(units.parsec)),
                                                                      N,
                                                                      r)

                    files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'

                    for f in files:
                        if 'sink' in f:
                            sink_particles = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)
                            #print sink_particles.time
                            #print sink_particles.x, sink_particles.y, sink_particles.z
                            #print '{0}/{1}'.format(filepath, f)
                            masses.append(numpy.mean(sink_particles.mass.value_in(units.MSun)))
            plots.append(numpy.mean(masses))

        pyplot.plot(Nsph, plots, lw=2, label='SFE={0}\%'.format(SFE[Mcloud.index(M)]))

    pyplot.xlabel(r'$N_\mathrm{SPH}$')
    pyplot.ylabel(r'$<M_\mathrm{sink}>$ [$M_{\odot}$]')
    pyplot.xticks(Nsph)
    pyplot.legend()
    pyplot.savefig('{0}/sink_mass.png'.format(save_path))
    pyplot.show()
    

def Nsinks_tff_vs_Nsph(path, save_path, Mcloud, Rcloud):

    Mcloud = [int(Mcloud.value_in(units.MSun))]#, 7500, 15000]
    SFE = [40, 25, 10]

    for M in Mcloud:

        rho_cloud = (M | units.MSun) / Rcloud ** 3
        tff = 1 / numpy.sqrt(constants.G * rho_cloud)
        print tff.value_in(units.Myr)

        Nsph = [4000, 8000, 16000, 32000]
        Nruns = [1, 2, 3, 4, 5, 6, 7, 8]

        plots = []

        for N in Nsph:
            masses = []
            for r in Nruns:
                if N * r <= 32000:
                    filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                                      M,
                                                                      int(Rcloud.value_in(units.parsec)),
                                                                      N,
                                                                      r)

                    files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'

                    for f in files:
                        if 'sink' in f:
                            sink_particles = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)
                            if sink_particles.get_timestamp() < tff:
                                print '{0}/{1}'.format(filepath, f)
                                print len(sink_particles)
            #plots.append(numpy.mean(masses))

        #pyplot.plot(Nsph, plots, lw=2, label='SFE={0}\%'.format(SFE[Mcloud.index(M)]))

    pyplot.xlabel(r'$N_\mathrm{SPH}$')
    pyplot.ylabel(r'$<M_\mathrm{sink}>$ [$M_{\odot}$]')
    pyplot.xticks(Nsph)
    pyplot.legend()
    pyplot.savefig('{0}/Nsinks.png'.format(save_path))
    #pyplot.show()


def mean_sink_size_vs_time(path, save_path, Mcloud, Rcloud):

    SFE = [40, 25, 10]

    Nsph = [4000, 8000, 16000, 32000]
    #Nruns = int(32000 / Nsph)

    for N in Nsph:
        Nruns = int(32000 / N)
        all_sizes = []
        all_times = []
        for r in range(1, Nruns + 1):
            filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                              int(Mcloud.value_in(units.MSun)),
                                                              int(Rcloud.value_in(units.parsec)),
                                                              N,
                                                              r)
            files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
            print filepath

            sizes = []
            times = []

            for f in files:
                if 'sink' in f:
                    sink_particles = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)

                    #print sink_particles.get_timestamp().value_in(units.Myr), \
                    #    numpy.mean(sink_particles.radius.value_in(units.parsec))

                    sizes.append(numpy.mean(sink_particles.radius.value_in(units.parsec)))
                    times.append(sink_particles.get_timestamp().value_in(units.Myr))

            sorted_times = numpy.sort(times)
            sorted_sizes = [x for _, x in sorted(zip(times, sizes))]

            #print sorted_times
            #print sorted_sizes

            all_times.append(sorted_times)
            all_sizes.append(sorted_sizes)

        pyplot.plot(numpy.mean(all_times, axis=0),
                    numpy.mean(all_sizes, axis=0),
                    label='N = {0}'.format(N))

    pyplot.xlabel(r'Time [Myr]')
    pyplot.ylabel(r'$<R_\mathrm{sink}>$ [pc]')
    pyplot.title('SFE = 40\%')
    pyplot.legend(loc="upper right")
    pyplot.savefig('{0}/mean_sink_size_vs_time.png'.format(save_path))
    pyplot.show()


def mean_sink_mass_vs_time(path, save_path, Mcloud, Rcloud):

    SFE = [40, 25, 10]

    Nsph = [4000, 8000, 16000, 32000]
    #Nruns = int(32000 / Nsph)

    for N in Nsph:
        Nruns = int(32000 / N)
        all_masses = []
        all_times = []
        for r in range(1, Nruns + 1):
            filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                              int(Mcloud.value_in(units.MSun)),
                                                              int(Rcloud.value_in(units.parsec)),
                                                              N,
                                                              r)
            files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
            print filepath

            masses = []
            times = []

            for f in files:
                if 'sink' in f:
                    sink_particles = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)

                    #print sink_particles.get_timestamp().value_in(units.Myr), \
                    #    numpy.mean(sink_particles.radius.value_in(units.parsec))

                    masses.append(numpy.mean(sink_particles.mass.value_in(units.MSun)))
                    times.append(sink_particles.get_timestamp().value_in(units.Myr))

            sorted_times = numpy.sort(times)
            sorted_masses = [x for _, x in sorted(zip(times, masses))]

            #print sorted_times
            #print sorted_sizes

            all_times.append(sorted_times)
            all_masses.append(sorted_masses)

        pyplot.plot(numpy.mean(all_times, axis=0),
                    numpy.mean(all_masses, axis=0),
                    label='N = {0}'.format(N))

    pyplot.xlabel(r'Time [Myr]')
    pyplot.ylabel(r'$<M_\mathrm{sink}>$ [$M_{\odot}$]')
    pyplot.title('SFE = 40\%')
    pyplot.legend(loc="lower right")
    pyplot.savefig('{0}/mean_sink_mass_vs_time.png'.format(save_path))
    pyplot.show()


def single_sink_mass_vs_time(path, save_path, Mcloud, Rcloud):

    from collections import defaultdict

    Nsph = [4000, 8000, 16000, 32000]

    for N in Nsph:

        Nruns = int(32000 / N)

        for r in range(1, Nruns + 1):
            filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                              int(Mcloud.value_in(units.MSun)),
                                                              int(Rcloud.value_in(units.parsec)),
                                                              N,
                                                              r)
            files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
            sink_files = [x for x in files if 'sink' in x]
            sink_files.sort(key=lambda f: int(filter(str.isdigit, f)))

            sink_masses = defaultdict(list)
            sink_times = defaultdict(list)

            for f in sink_files:
                sinks = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)
                print sinks.tff.value_in(units.Myr)

                for s in sinks:
                    sink_masses[s.key].append(s.mass.value_in(units.MSun))
                    sink_times[s.key].append(sinks.get_timestamp().value_in(units.Myr))

            for key, val in sink_masses.items():
                if len(val) == 1:
                    pyplot.plot(sink_times[key], val, ls='-', marker='o', markersize=4)
                else:
                    pyplot.plot(sink_times[key], val, ls='-')

            pyplot.xlabel(r'Time [Myr]')
            pyplot.ylabel(r'Sink mass [$M_{\odot}$]')
            pyplot.title(r'SFE = 40\%, $N_{{SPH}} = {0}, M_{{cloud}} = {1} M_{{\odot}}$'.format(N,
                                                                                                int(Mcloud.value_in(units.MSun))))
            figname = '{0}/M{1}MSun_R{2}pc_N{3}_r{4}.png'.format(save_path,
                                                                 int(Mcloud.value_in(units.MSun)),
                                                                 int(Rcloud.value_in(units.parsec)),
                                                                 N,
                                                                 r)
            pyplot.savefig(figname)
            pyplot.show()


def total_sink_mass_vs_time(path, save_path, Mcloud, Rcloud):
    Nsph = [4000, 8000, 16000, 32000]
    #Nruns = int(32000 / Nsph)

    for N in Nsph:
        Nruns = int(32000 / N)
        all_masses = []
        all_times = []
        for r in range(1, Nruns + 1):
            filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                              int(Mcloud.value_in(units.MSun)),
                                                              int(Rcloud.value_in(units.parsec)),
                                                              N,
                                                              r)
            files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
            print filepath

            masses = []
            times = []

            for f in files:
                if 'sink' in f:
                    sink_particles = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)

                    #print sink_particles.get_timestamp().value_in(units.Myr), \
                    #    numpy.mean(sink_particles.radius.value_in(units.parsec))

                    masses.append(sink_particles.mass.sum().value_in(units.MSun))
                    times.append(sink_particles.get_timestamp().value_in(units.Myr))

            sorted_times = numpy.sort(times)
            sorted_masses = [x for _, x in sorted(zip(times, masses))]

            #print sorted_times
            #print sorted_sizes

            all_times.append(sorted_times)
            all_masses.append(sorted_masses)

        pyplot.plot(numpy.mean(all_times, axis=0),
                    numpy.mean(all_masses, axis=0),
                    label='N = {0}'.format(N))

    pyplot.xlabel(r'Time [Myr]')
    pyplot.ylabel(r'$\sum M_\mathrm{sink}$ [$M_{\odot}$]')
    pyplot.title('SFE = 40\%')
    pyplot.legend(loc="lower right")
    pyplot.savefig('{0}/total_sink_mass_vs_time.png'.format(save_path))
    pyplot.show()


def Nsinks_vs_time(path, save_path, Mcloud, Rcloud):

    SFE = [40, 25, 10]

    Nsph = [4000, 8000, 16000, 32000]
    #Nruns = int(32000 / Nsph)

    for N in Nsph:
        Nruns = int(32000 / N)
        all_sinks = []
        all_times = []
        for r in range(1, Nruns + 1):
            filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                              int(Mcloud.value_in(units.MSun)),
                                                              int(Rcloud.value_in(units.parsec)),
                                                              N,
                                                              r)
            files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
            print filepath

            sinks = []
            times = []

            for f in files:
                if 'sink' in f:
                    sink_particles = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)
                    #print sink_particles.get_timestamp().value_in(units.Myr), \
                    #    numpy.mean(sink_particles.radius.value_in(units.parsec))

                    sinks.append(len(sink_particles))
                    times.append(sink_particles.get_timestamp().value_in(units.Myr))

            sorted_times = numpy.sort(times)
            sorted_sinks = [x for _, x in sorted(zip(times, sinks))]

            #print sorted_times
            #print sorted_sizes

            all_times.append(sorted_times)
            all_sinks.append(sorted_sinks)

        pyplot.plot(numpy.mean(all_times, axis=0),
                    numpy.mean(all_sinks, axis=0),
                    label='N = {0}'.format(N))

    pyplot.xlabel(r'Time [Myr]')
    pyplot.ylabel(r'$N_\mathrm{sink}$')
    pyplot.title('SFE = 40\%')
    pyplot.legend(loc="best")
    pyplot.savefig('{0}/Nsinks_vs_time.png'.format(save_path))
    pyplot.show()


def sink_location_vs_time(path, save_path, Mcloud, Rcloud):

    SFE = [40, 25, 10]

    Nsph = [4000, 8000, 16000, 32000]
    #Nruns = int(32000 / Nsph)

    cmap = pyplot.cm.get_cmap('Oranges', 100)  # PiYG
    colors = []

    for i in range(cmap.N):
        rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
        colors.append(pyplot.cm.colors.rgb2hex(rgb))

    fig, axes = pyplot.subplots(2, 2, figsize=(14, 14))

    # Plot "cloud"
    cloud = pyplot.Circle((0, 0), 1.0, color='k', alpha=0.2, fill=False)
    #for ax in axes.flatten():
    #    ax.set_aspect('equal')
    #    ax.add_artist(cloud)

    # Subplots locations
    subpl = {4000: (0, 0), 8000: (0, 1), 16000: (1, 0), 32000: (1, 1)}

    for N in Nsph:
        filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                          int(Mcloud.value_in(units.MSun)),
                                                          int(Rcloud.value_in(units.parsec)),
                                                          N,
                                                          1)
        files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
        print filepath

        locs_x = []
        locs_y = []
        times = []

        for f in files:
            if 'sink' in f:
                sink_particles = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)

                locs_x.append(numpy.mean(sink_particles.x.value_in(units.parsec)))
                locs_y.append(numpy.mean(sink_particles.y.value_in(units.parsec)))
                times.append(sink_particles.get_timestamp().value_in(units.Myr))
                #pyplot.plot(sink_particles.get_timestamp().value_in(units.Myr),
                #            numpy.mean(sink_particles.x.value_in(units.parsec)),
                #            'ro')

        sorted_times = numpy.sort(times)
        sorted_x = [x for _, x in sorted(zip(times, locs_x))]
        sorted_y = [x for _, x in sorted(zip(times, locs_y))]

        print len(sorted_times)

        axes[subpl[N][0], subpl[N][1]].scatter(sorted_x, sorted_y, c=colors)
        axes[subpl[N][0], subpl[N][1]].set_aspect('equal')
        #axes[subpl[N][0], subpl[N][1]].add_artist(cloud)
        axes[subpl[N][0], subpl[N][1]].set_xlim([-1.0, 1.0])
        axes[subpl[N][0], subpl[N][1]].set_ylim([-1.0, 1.0])
        axes[subpl[N][0], subpl[N][1]].set_xlabel(r'x [pc]')
        axes[subpl[N][0], subpl[N][1]].set_ylabel(r'y [pc]')

    # Had to do this because pyplot was being weird about the colorbar's cmap
    norm = pyplot.Normalize(sorted_times[0], sorted_times[-1])
    sm = pyplot.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    fig.subplots_adjust(right=0.8, hspace=0.5)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(sm, cax=cbar_ax, label='Time [Myr]')#, fraction=0.046, pad=0.04)

    pyplot.suptitle('SFE = 40\%')
    pyplot.savefig('{0}/sinks_location.png'.format(save_path))
    pyplot.show()


def stars_locations(path, save_path, Rcloud, Nsph, Mcloud):
    filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                      int(Mcloud.value_in(units.MSun)),
                                                      int(Rcloud.value_in(units.parsec)),
                                                      Nsph,
                                                      1)
    filepath = path
    files = os.listdir(filepath)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
    stars_files = [x for x in files if 'stars' in x]
    stars_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    #stars = read_set_from_file('{0}/{1}'.format(filepath, stars_files[-1]), "hdf5", close_file=True)

    i = 0

    for sf in stars_files:
        fig = pyplot.figure(figsize=(8, 8))
        stars = read_set_from_file('{0}/{1}'.format(filepath, sf), "hdf5", close_file=True)
        print stars_files[-1]
        print len(stars)
        pyplot.scatter(stars.x.value_in(units.parsec),
                       stars.y.value_in(units.parsec))

        ax = fig.gca()
        ax.set_aspect('equal')
        ax.set_xlim([-1.0, 1.0])
        ax.set_ylim([-1.0, 1.0])
        ax.set_xlabel(r'x [pc]')
        ax.set_ylabel(r'y [pc]')

        pyplot.savefig('{0}/{1}.png'.format(save_path, i))
        i += 1


def star_formation_movie(path, save_path, Rcloud, Nsph, Mcloud):
    filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                      int(Mcloud.value_in(units.MSun)),
                                                      int(Rcloud.value_in(units.parsec)),
                                                      Nsph,
                                                      1)
    filepath = path
    files = os.listdir(filepath)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'

    sink_files = [x for x in files if 'sink' in x]
    sink_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    star_files = [x for x in files if 'stars' in x]
    star_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    #stars = read_set_from_file('{0}/{1}'.format(filepath, stars_files[-1]), "hdf5", close_file=True)

    i = 0

    prev_stars = None
    times = []
    Nstars = []
    Nstars_hm = []

    for s in sink_files:
        print s
        sink_file_number = s.split('_')[-1].split('.')[0][1:]

        fig = pyplot.figure(figsize=(8, 8))
        sinks = read_set_from_file('{0}/{1}'.format(filepath, s), "hdf5", close_file=True)
        pyplot.scatter(sinks.x.value_in(units.parsec),
                       sinks.y.value_in(units.parsec),
                       alpha=0.5)

        for sf in star_files:
            star_file_number = sf.split('_')[-1].split('.')[0][1:]

            if star_file_number == sink_file_number:
                stars = read_set_from_file('{0}/{1}'.format(filepath, sf), "hdf5", close_file=True)
                if prev_stars is None:
                    prev_stars = stars.key
                else:
                    # Is there a way to do this??
                    #old_stars = stars[stars.key in prev_stars.key]
                    #new_stars = stars[stars.key not in prev_stars.key]
                    #print old_stars
                    #prev_stars = stars.key
                    #print "old stars: {0}, new stars: {1}".format(len(old_stars), len(new_stars))
                    for star in stars:
                        if star.key in prev_stars:  # old star
                            pyplot.scatter(star.x.value_in(units.parsec),
                                           star.y.value_in(units.parsec), marker="*", color='blue')
                        else:  # new star
                            pyplot.scatter(star.x.value_in(units.parsec),
                                           star.y.value_in(units.parsec), marker="*", color='black')
                    prev_stars = stars.key
                times.append(stars.get_timestamp().value_in(units.Myr))
                Nstars.append(len(stars))
                Nstars_hm.append(len(stars[stars.mass >= 1.9 | units.MSun]))

        ax = fig.gca()
        ax.set_aspect('equal')
        ax.set_xlim([-2.0, 2.0])
        ax.set_ylim([-2.0, 2.0])
        ax.set_xlabel(r'x [pc]')
        ax.set_ylabel(r'y [pc]')

        pyplot.savefig('{0}/{1}.png'.format(save_path, i))
        i += 1

    fig = pyplot.figure(figsize=(8, 8))
    pyplot.plot(times, Nstars, label="All stars")
    pyplot.plot(times, Nstars_hm, label=r"$M_* \geq 1.9 M_{\odot}$")
    pyplot.legend(loc='lower right')

    ax = fig.gca()
    #ax.set_aspect('equal')
    #ax.set_xlim([-1.0, 1.0])
    #ax.set_ylim([-1.0, 1.0])
    ax.set_xlabel(r'Time [Myr]')
    ax.set_ylabel(r'$N_*$')

    pyplot.savefig('{0}/Nstars.png'.format(save_path))


def main(path, save_path, tend, dt_diag, Ncloud, Mcloud, Rcloud):
    # My own style sheet, comment out if not needed
    pyplot.style.use('paper')

    #stars_locations(path, save_path, Rcloud, 4000, Mcloud)
    #star_formation_movie(path, save_path, Rcloud, 4000, Mcloud)

    #mean_sink_size_vs_Nsph(path, save_path, Mcloud, Rcloud)
    #mean_sink_mass_vs_Nsph(path, save_path, Mcloud, Rcloud)
    #Nsinks_tff_vs_Nsph(path, save_path, Mcloud, Rcloud)

    #Nsinks_vs_time(path, save_path, Mcloud, Rcloud)
    #mean_sink_size_vs_time(path, save_path, Mcloud, Rcloud)
    #mean_sink_mass_vs_time(path, save_path, Mcloud, Rcloud)
    #total_sink_mass_vs_time(path, save_path, Mcloud, Rcloud)
    #sink_location_vs_time(path, save_path, Mcloud, Rcloud)

    single_sink_mass_vs_time(path, save_path, Mcloud, Rcloud)



def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-p", dest="path",
                      default=".",
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

    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    numpy.random.seed(3141)
    main(**o.__dict__)

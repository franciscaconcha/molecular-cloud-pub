import numpy
from matplotlib import pyplot
import os

from amuse.lab import *


def Nsph_vs_mean_sink_size(path, save_path, Rcloud):

    Mcloud = [4000, 7000, 15000]
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
                            #print sink_particles.time
                            #print sink_particles.x, sink_particles.y, sink_particles.z
                            print '{0}/{1}'.format(filepath, f)
                            sizes.append(numpy.mean(sink_particles.radius.value_in(units.parsec)))
            plots.append(numpy.mean(sizes))

        pyplot.plot(Nsph, plots, lw=2, label='SFE={0}\%'.format(SFE[Mcloud.index(M)]))

    pyplot.xlabel(r'$N_\mathrm{SPH}$')
    pyplot.ylabel(r'$<R_\mathrm{sink}>$ [pc]')
    pyplot.xticks(Nsph)
    pyplot.legend()
    pyplot.savefig('{0}/sink_size.png'.format(save_path))
    pyplot.show()


def Nsph_vs_mean_sink_mass(path, save_path, Rcloud):

    Mcloud = [4000, 7000, 15000]
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
                            print '{0}/{1}'.format(filepath, f)
                            masses.append(numpy.mean(sink_particles.mass.value_in(units.MSun)))
            plots.append(numpy.mean(masses))

        pyplot.plot(Nsph, plots, lw=2, label='SFE={0}\%'.format(SFE[Mcloud.index(M)]))

    pyplot.xlabel(r'$N_\mathrm{SPH}$')
    pyplot.ylabel(r'$<M_\mathrm{sink}>$ [$M_{\odot}$]')
    pyplot.xticks(Nsph)
    pyplot.legend()
    pyplot.savefig('{0}/sink_mass.png'.format(save_path))
    pyplot.show()


def time_vs_mean_sink_size(path, save_path, Rcloud):

    Mcloud = 4000#, 7500, 15000]
    SFE = [40, 25, 10]

    Nsph = [4000, 8000, 16000, 32000]
    #Nruns = int(32000 / Nsph)

    for N in Nsph:
        Nruns = int(32000 / N)
        all_sizes = []
        all_times = []
        for r in range(1, Nruns + 1):
            filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                              Mcloud,
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
    pyplot.legend()
    pyplot.savefig('{0}/time_vs_mean_sink_size.png'.format(save_path))
    pyplot.show()


def time_vs_mean_sink_mass(path, save_path, Rcloud):

    Mcloud = 7500#, 15000]
    SFE = [40, 25, 10]

    Nsph = [4000, 8000, 16000, 32000]
    #Nruns = int(32000 / Nsph)

    for N in Nsph:
        Nruns = int(32000 / N)
        all_masses = []
        all_times = []
        for r in range(1, Nruns + 1):
            filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                              Mcloud,
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
    pyplot.legend()
    pyplot.savefig('{0}/time_vs_mean_sink_mass.png'.format(save_path))
    pyplot.show()


def time_vs_sink_location(path, save_path, Rcloud):

    Mcloud = 4000#, 15000]
    SFE = [40, 25, 10]

    Nsph = [4000]#[4000, 8000, 16000, 32000]
    #Nruns = int(32000 / Nsph)

    fig = pyplot.figure()

    cmap = pyplot.cm.get_cmap('Oranges', 100)  # PiYG
    colors = []

    for i in range(cmap.N):
        rgb = cmap(i)[:100]  # will return rgba, we take only first 3 so we get rgb
        colors.append(pyplot.cm.colors.rgb2hex(rgb))
    i = 0  # For colors

    for N in Nsph:
        Nruns = int(32000 / N)
        all_locs_x = []
        all_locs_y = []
        all_times = []
        Nruns = 2
        for r in range(1, Nruns + 1):
            filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                              Mcloud,
                                                              int(Rcloud.value_in(units.parsec)),
                                                              N,
                                                              r)
            files = os.listdir(filepath) #= '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'
            print filepath

            locs_x = []
            locs_y = []
            times = []

            for f in files:
                if 'sink' in f:
                    sink_particles = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)

                    #print sink_particles.get_timestamp().value_in(units.Myr), \
                    #    numpy.mean(sink_particles.radius.value_in(units.parsec))

                    locs_x.append(numpy.mean(sink_particles.x.value_in(units.parsec)))
                    locs_y.append(numpy.mean(sink_particles.y.value_in(units.parsec)))
                    times.append(sink_particles.get_timestamp().value_in(units.Myr))

            sorted_times = numpy.sort(times)
            sorted_x = [x for _, x in sorted(zip(times, locs_x))]
            sorted_y = [x for _, x in sorted(zip(times, locs_y))]

            #print sorted_times
            #print sorted_sizes

            #all_times.append(sorted_times)
            #all_masses.append(sorted_masses)

            print len(sorted_times)

            all_times.append(sorted_times)
            all_locs_x.append(sorted_x)
            all_locs_y.append(sorted_y)

            #for t in sorted_times:
            #    s = pyplot.scatter(sorted_x[list(sorted_times).index(t)],
            #                       sorted_y[list(sorted_times).index(t)],
            #                       color=colors[i], label=t, cmap=cmap)#,
            #    i += 1
                    #label='N = {0}'.format(N))

    all_locs_x = numpy.hstack(all_locs_x)
    all_locs_y = numpy.hstack(all_locs_y)

    norm = pyplot.Normalize(sorted_times[0], sorted_times[-1])
    s = pyplot.scatter(all_locs_x, all_locs_y, c=colors)

    sm = pyplot.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm)

    pyplot.xlabel(r'x [pc]')
    pyplot.ylabel(r'y [pc]')
    pyplot.title('SFE = 40\%')
    #cb = pyplot.colorbar(s)
    cbar.set_label('Time [Myr]')
    pyplot.xlim([-1.0, 1.0])
    pyplot.ylim([-1.0, 1.0])
    #pyplot.savefig('{0}/time_vs_mean_sink_mass.png'.format(save_path))
    pyplot.show()


def main(path, save_path, tend, dt_diag, Ncloud, Mcloud, Rcloud):
    import os
    pyplot.style.use('paper')

    #Nsph_vs_mean_sink_size(path, save_path, Rcloud)
    #Nsph_vs_mean_sink_mass(path, save_path, Rcloud)

    #time_vs_mean_sink_size(path, save_path, Rcloud)
    #time_vs_mean_sink_mass(path, save_path, Rcloud)
    time_vs_sink_location(path, save_path, Rcloud)

    #sink_motion(path, save_path, Rcloud)


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

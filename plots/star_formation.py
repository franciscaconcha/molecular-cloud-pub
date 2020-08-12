import numpy
from matplotlib import pyplot
import os

from amuse.lab import *


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

    # stars = read_set_from_file('{0}/{1}'.format(filepath, stars_files[-1]), "hdf5", close_file=True)

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
    """filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                      int(Mcloud.value_in(units.MSun)),
                                                      int(Rcloud.value_in(units.parsec)),
                                                      Nsph,
                                                      1)"""
    filepath = path
    files = os.listdir(filepath)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'

    sink_files = [x for x in files if 'sink_particles' in x]
    sink_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    # star_files = [x for x in files if 'stars_particles' in x]
    # star_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    stars_file = [x for x in files if 'gravity_stars' in x]
    stars = read_set_from_file('{0}/{1}'.format(filepath, stars_file[0]), "hdf5", close_file=True)
    print len(stars)

    i = 0

    prev_stars = None
    times = []
    Nstars = []
    Nstars_hm = []

    tprev = 0.0 | units.Myr

    # Find first timestamp where stars are formed
    tmin = 20.0 | units.Myr
    for s in stars:
        if s.tborn < tmin:
            tmin = s.tborn

    for s in sink_files:
        print s
        sink_file_number = s.split('_')[-1].split('.')[0][1:]

        fig = pyplot.figure(figsize=(8, 8))
        sinks = read_set_from_file('{0}/{1}'.format(filepath, s), "hdf5", close_file=True)
        pyplot.scatter(sinks.x.value_in(units.parsec),
                       sinks.y.value_in(units.parsec),
                       alpha=0.5)
        time = sinks.get_timestamp()

        if time >= tmin:  # Stars are forming now
            # All the stars formed before tprev
            old_stars = stars[stars.tborn <= tprev]
            pyplot.scatter(old_stars.x.value_in(units.parsec),
                           old_stars.y.value_in(units.parsec), marker="*", color='blue')

            # Finding stars formed between tprev and time, these are the new stars
            prev_stars = stars[stars.tborn > tprev]
            new_stars = prev_stars[prev_stars.tborn <= time]
            pyplot.scatter(new_stars.x.value_in(units.parsec),
                           new_stars.y.value_in(units.parsec), marker="*", color='black')
            tprev = time

        times.append(time.value_in(units.Myr))
        Nstars.append(len(old_stars) + len(new_stars))
        Nstars_hm.append(len(old_stars[old_stars.stellar_mass >= 1.9 | units.MSun]) +
                         len(new_stars[new_stars.stellar_mass >= 1.9 | units.MSun]))

        """for sf in stars_files:
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
                Nstars_hm.append(len(stars[stars.stellar_mass >= 1.9 | units.MSun]))"""

        ax = fig.gca()
        ax.set_aspect('equal')
        ax.set_xlim([-2.0, 2.0])
        ax.set_ylim([-2.0, 2.0])
        ax.set_xlabel(r'x [pc]')
        ax.set_ylabel(r'y [pc]')

        pyplot.savefig('{0}/{1}.png'.format(save_path, i))
        pyplot.close()
        i += 1

    fig = pyplot.figure(figsize=(8, 8))
    pyplot.plot(times, Nstars, label="All stars")
    pyplot.plot(times, Nstars_hm, label=r"$M_* \geq 1.9 M_{\odot}$")
    pyplot.legend(loc='lower right')

    ax = fig.gca()
    # ax.set_aspect('equal')
    # ax.set_xlim([-1.0, 1.0])
    # ax.set_ylim([-1.0, 1.0])
    ax.set_xlabel(r'Time [Myr]')
    ax.set_ylabel(r'$N_*$')

    pyplot.savefig('{0}/Nstars.png'.format(save_path))


def a_vs_e(path, save_path, N):
    filepath = path
    files = os.listdir(filepath)  # = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'

    sink_files = [x for x in files if 'sink_particles' in x]
    sink_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    # star_files = [x for x in files if 'stars_particles' in x]
    # star_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    stars_file = [x for x in files if 'gravity_stars' in x]
    stars = read_set_from_file('{0}/{1}'.format(filepath, stars_file[0]), "hdf5", close_file=True)
    print len(stars)

    stars.mass = stars.stellar_mass

    converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(),
                                         stars[1].position.length())

    binaries = stars.get_binaries()
    print len(binaries)

    for b in binaries:
        kep = Kepler(converter)
        kep.initialize_from_particles(b)
        a, e = kep.get_elements()

        pyplot.scatter(a, e)
    pyplot.show()


def final_imf(path, save_path, Mcloud, Rcloud, N, r=1):
    fig = pyplot.figure()
    ax = pyplot.gca()
    filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                      int(Mcloud.value_in(units.MSun)),
                                                      int(Rcloud.value_in(units.parsec)),
                                                      N,
                                                      r)
    files = os.listdir(filepath)
    star_files = [x for x in files if 'star' in x]
    star_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    final_stars = read_set_from_file('{0}/{1}'.format(filepath, star_files[-1]), "hdf5", close_file=True)
    star_masses = final_stars.mass.value_in(units.MSun)

    # Create a Kroupa 2001 distribution with the same number of stars, to compare
    IMF_masses = new_kroupa_mass_distribution(len(final_stars),
                                              mass_max=max(star_masses) | units.MSun).value_in(units.MSun)

    kroupa_limits = [0.01, 0.08, 0.5, 1.0, 1.9, max(star_masses)]  # Added the 1.9 MSun 'bin' for photoevap limit

    ax.hist(IMF_masses, bins=kroupa_limits, label='Kroupa IMF', histtype=u'step', edgecolor='r', lw=3)
    ax.hist(star_masses, bins=kroupa_limits, label='Simulation', histtype=u'step', edgecolor='k', lw=3)

    # 1.9 MSun limit, for photoevaporation
    ax.axvline(1.9,
               lw=3,
               ls='--',
               c='blue')
    ax.text(0.666,
            0.5,
            r'$N_*(M_* >= 1.9 M_{{\odot}}$) = {0}$'.format(len(star_masses[star_masses >= 1.9])),
            color='blue',
            transform=ax.transAxes)

    ax.set_xlabel(r'$M_*$ [$\mathrm{M}_{\odot}$]')
    ax.set_ylabel(r'$N_*$')
    ax.set_title(r'Total $N_*$ = {0}'.format(len(final_stars)))
    ax.set_xscale('log')
    pyplot.legend(loc='best')
    pyplot.show()
    pyplot.savefig('{0}/IMF_vs_simulation.png'.format(save_path))


def Nstars_vs_time(path, save_path, Mcloud, Rcloud, N, r=1):
    filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                      int(Mcloud.value_in(units.MSun)),
                                                      int(Rcloud.value_in(units.parsec)),
                                                      N,
                                                      r)
    files = os.listdir(filepath)

    star_files = [x for x in files if 'stars' in x]
    star_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    Nstars, times = [], []

    for f in star_files:
        stars = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)
        time = stars.get_timestamp().value_in(units.Myr)

        times.append(time)
        Nstars.append(len(stars))

    import decimal
    # Need to do this to round down the last time stamp properly
    d = float(decimal.Decimal(times[0]).quantize(decimal.Decimal('.1'), rounding=decimal.ROUND_DOWN))

    earlier_times = numpy.arange(0., d, 0.1)
    times = numpy.concatenate((earlier_times, times))

    Nstars = numpy.pad(Nstars, (len(times) - len(Nstars), 0), 'constant')

    pyplot.plot(times, Nstars, lw=3)

    pyplot.xlabel('Time [Myr]')
    pyplot.ylabel(r'$N_*$')
    pyplot.title(r'Final $N_*$ = {0}'.format(len(stars)))  # Current stars is the last file
    pyplot.savefig('{0}/Nstars_vs_time.png'.format(path))
    pyplot.show()


def Mstars_vs_time(path, save_path, Mcloud, Rcloud, N, r=1):
    filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                      int(Mcloud.value_in(units.MSun)),
                                                      int(Rcloud.value_in(units.parsec)),
                                                      N,
                                                      r)
    files = os.listdir(filepath)

    star_files = [x for x in files if 'stars' in x]
    star_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    for f in star_files:
        stars = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)
        star_masses = stars.mass.value_in(units.MSun)
        time = stars.get_timestamp().value_in(units.Myr)

        pyplot.scatter(time * numpy.ones(star_masses.shape), star_masses,
                       marker='*',
                       edgecolors='k',
                       facecolors='k',
                       alpha=0.5)

    pyplot.axhline(1.9, color='r')
    pyplot.text(0.01, 2, r'$1.9 M_{\odot}$', color='r')

    pyplot.xlim([0.0, 0.8])
    pyplot.xlabel('Time [Myr]')
    pyplot.ylabel(r'$M_*$')
    pyplot.title(r'Final $N_*$ = {0}'.format(len(stars)))  # Current stars is the last file
    pyplot.savefig('{0}/Mstars_vs_time.png'.format(path))
    pyplot.show()


def stars_vs_time(path, save_path, Mcloud, Rcloud, N, r=1):
    fig, axes = pyplot.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)

    filepath = '{0}/M{1}MSun_R{2}pc_N{3}/{4}/'.format(path,
                                                      int(Mcloud.value_in(units.MSun)),
                                                      int(Rcloud.value_in(units.parsec)),
                                                      N,
                                                      r)
    files = os.listdir(filepath)

    star_files = [x for x in files if 'stars' in x]
    star_files.sort(key=lambda f: int(filter(str.isdigit, f)))

    Nstars, times = [], []

    for f in star_files:
        stars = read_set_from_file('{0}/{1}'.format(filepath, f), "hdf5", close_file=True)
        star_masses = stars.mass.value_in(units.MSun)
        time = stars.get_timestamp().value_in(units.Myr)

        times.append(time)
        Nstars.append(len(stars))

        axes[1].scatter(time * numpy.ones(star_masses.shape), star_masses,
                        marker='*',
                        edgecolors='k',
                        facecolors='k',
                        alpha=0.5)

    axes[1].axhline(1.9, color='r', ls='--', alpha=0.5)
    axes[1].text(0.01, 2, r'$1.9 M_{\odot}$', color='r', alpha=0.5)

    import decimal
    # Need to do this to round down the last time stamp properly
    d = float(decimal.Decimal(times[0]).quantize(decimal.Decimal('.01'), rounding=decimal.ROUND_DOWN))

    dt = times[0] - d  # Have to do this so that earlier_times[-1] is always earlier than times[0]
    earlier_times = numpy.arange(0., d + dt, dt)
    times = numpy.concatenate((earlier_times, times))

    Nstars = numpy.pad(Nstars, (len(times) - len(Nstars), 0), 'constant')

    axes[0].plot(times, Nstars, lw=3, c='#348ABD')

    axes[0].set_xlabel('Time [Myr]')
    axes[0].set_ylabel(r'$N_*$')
    axes[0].set_title(r'Final $N_*$ = {0}'.format(len(stars)))  # Current stars is the last file

    axes[1].set_xlim([0.0, 0.8])
    axes[1].set_xlabel('Time [Myr]')
    axes[1].set_ylabel(r'$M_*$')
    pyplot.savefig('{0}/stars_vs_time.png'.format(path))
    pyplot.show()


def main(path, save_path, tend, dt_diag, Ncloud, Mcloud, Rcloud):
    # My own style sheet, comment out if not needed
    pyplot.style.use('paper')
    pyplot.close('all')

    # stars_locations(path, save_path, Rcloud, 4000, Mcloud)
    star_formation_movie(path, save_path, Rcloud, 24000, Mcloud)
    #a_vs_e(path, save_path, Ncloud)

    # final_imf(path, save_path, Mcloud, Rcloud, Ncloud)

    # Nstars_vs_time(path, save_path, Mcloud, Rcloud, Ncloud)
    # Mstars_vs_time(path, save_path, Mcloud, Rcloud, Ncloud)
    # stars_vs_time(path, save_path, Mcloud, Rcloud, Ncloud)


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
                      default=4000,
                      type="float",
                      help="number of gas particles.")
    result.add_option("--Mcloud", dest="Mcloud",
                      unit=units.MSun,
                      type="float",
                      default=10000 | units.MSun,
                      help="cloud mass")
    result.add_option("--Rcloud", dest="Rcloud",
                      unit=units.parsec,
                      type="float",
                      default=2 | units.parsec,
                      help="cloud size")

    return result


if __name__ in ("__main__", "__plot__"):
    o, arguments = new_option_parser().parse_args()
    numpy.random.seed(3141)
    main(**o.__dict__)

import numpy
import scipy.interpolate


# Set to true if your environment can run amuse
AMUSE_ENABLED = True

if AMUSE_ENABLED:
    from amuse.lab import *
    G0 = 1.6e-3 * units.erg / units.s / units.cm**2.

class FRIED_interpolator:
    """ FRIED interpolator object.
    Performs linear interpolation on the FRIED grid.
    As there is a degeneracy between disk radius, disk mass and disk outer density, we neglect disk outer density.
    Interpolation can optionally be done in lin or log space, separately for each axis
    """
    def __init__(self,
                 folder='data',
                 logMstar=True,
                 logF=True,
                 logMdisk=True,
                 logRdisk=True,
                 verbosity=False):
        """
        Initialize FRIED_interpolator object

        :param folder: directory of FRIED grid data file [string]
        :param logMstar: construct grid of the log of the host star mass (else linear) [boolean]
        :param logF: construct grid of the log of the incident FUV field (else linear) [boolean]
        :param logMdisk: construct grid of the log of the disk mass (else linear) [boolean]
        :param logRdisk: construct grid of the log of the disk radius (else linear) [boolean]
        """

        Mstar_grid, F_grid, Mdisk_grid, Rdisk_grid = numpy.loadtxt(folder + '/friedgrid.dat',
                                                                   usecols=(0, 1, 2, 4),
                                                                   unpack=True)

        self._logMdot_grid = numpy.loadtxt(folder + '/friedgrid.dat', usecols=(5,))

        self._logMstar = logMstar
        self._logF = logF
        self._logMdisk = logMdisk
        self._logRdisk = logRdisk

        if self._logMstar:
            Mstar_grid = numpy.log10(Mstar_grid)
        if self._logF:
            F_grid = numpy.log10(F_grid)
        if self._logMdisk:
            Mdisk_grid = numpy.log10(Mdisk_grid)
        if self._logRdisk:
            Rdisk_grid = numpy.log10(Rdisk_grid)

        self._grid = numpy.array([Mstar_grid,
                                  F_grid,
                                  Mdisk_grid,
                                  Rdisk_grid,
                                  ]).T

        self._interpolator = scipy.interpolate.LinearNDInterpolator(self._grid, self._logMdot_grid)
        self._backup_interpolator = scipy.interpolate.NearestNDInterpolator(self._grid, self._logMdot_grid)

        self.verbosity = verbosity

    def interp(self, Mstar, F, Mdisk, Rdisk):
        """
        Compute mass loss rate at N positions on the grid
        Input and output is in standard FRIED units
        Lin/log space interpolation is handled in code

        :param Mstar: host star mass [1D float array, shape (N)]
        :param F: host star mass [1D float array, shape (N)]
        :param Mdisk: disk mass [1D float array, shape (N)]
        :param Rdisk: disk radius [1D float array, shape (N)]

        :return: mass loss rate for each set of parameters [1D float array, shape (N)]
        """

        if self._logMstar:
            Mstar = numpy.log10(Mstar)
        if self._logF:
            F = numpy.log10(F)
        if self._logMdisk:
            Mdisk = numpy.log10(Mdisk)
        if self._logRdisk:
            Rdisk = numpy.log10(Rdisk)

        x_i = numpy.array([Mstar,
                           F,
                           Mdisk,
                           Rdisk
                           ]).T

        logMdot = self._interpolator(x_i)

        if numpy.isnan(logMdot):
            logMdot = self._backup_interpolator(x_i)

        return 10.**logMdot

    if AMUSE_ENABLED:
        def interp_amuse(self, Mstar, F, Mdisk, Rdisk):
            """
            Compute mass loss rate at N positions on the grid
            Input and output is in amuse units
            Lin/log space interpolation is handled in code

            :param Mstar: host star mass [1D mass unit scalar array, shape (N)]
            :param F: host star mass [1D flux unit scalar array, shape (N)]
            :param Mdisk: disk mass [1D mass unit scalar array, shape (N)]
            :param Rdisk: disk radius [1D length unit scalar array, shape (N)]

            :raturn: mass loss rate for each set of parameters [1D mass flux per time unit scalar array, shape (N)]
            """

            if self._logMstar:
                Mstar = numpy.log10(Mstar.value_in(units.MSun))
            else:
                Mstar = Mstar.value_in(units.MSun)
            if self._logF:
                F = numpy.log10(F.value_in(G0))
            else:
                F = F.value_in(G0)
            if self._logMdisk:
                Mdisk = numpy.log10(Mdisk.value_in(units.MJupiter))
            else:
                Mdisk = Mdisk.value_in(units.MJupiter)
            if self._logRdisk:
                Rdisk = numpy.log10(Rdisk.value_in(units.AU))
            else:
                Rdisk = Rdisk.value_in(units.AU)

            x_i = numpy.array([Mstar,
                               F,
                               Mdisk,
                               Rdisk
                               ]).T

            logMdot = self._interpolator(x_i)

            if numpy.isnan(logMdot):
                logMdot = self._backup_interpolator(x_i)
                if self.verbosity:
                    print ("[WARNING] Point outside interpolation domain, falling back to nearest neighbour")

            return 10.**logMdot | units.MSun / units.yr


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import time

    N = 1000

    Mdisk = numpy.logspace(-3., 3., num=N)
    Mstar = 1.*numpy.ones(len(Mdisk))
    F = 1000.*numpy.ones(len(Mdisk))
    Rdisk = 100.*numpy.ones(len(Mdisk))


    print ('Interpolations in log space:')
    start = time.time()

    interpolator = FRIED_interpolator()

    end = time.time()

    print ('Building the interpolator for a {a} point grid took {b} s'.format(a=interpolator._grid.shape[0], b=end-start))


    start = time.time()

    Mdot = interpolator.interp(Mstar, F, Mdisk, Rdisk)

    end = time.time()

    print ('Performing {a} interpolations took {b} s, \nor {c} s per interpolation'.format(a=N, b=end-start, c=(end-start)/N))


    on_grid_mask = (interpolator._grid[:,0] == numpy.log10(Mstar[0]))*\
                   (interpolator._grid[:,1] == numpy.log10(F[0]))*\
                   (interpolator._grid[:,3] == numpy.log10(Rdisk[0]))

    Mdisk_on_grid = interpolator._grid[ on_grid_mask, 2 ]

    plt.figure()
    plt.scatter(Mdisk_on_grid, interpolator._logMdot_grid[ on_grid_mask ])
    plt.plot(numpy.log10(Mdisk), numpy.log10(Mdot))

    plt.xlabel('log10 Disk Mass (MJup)')
    plt.ylabel('log10 Disk Mass Loss Rate (MSun/yr)')

    plt.title('Interpolations along disk mass slice, in log space')

    print ('Interpolations in lin space:')
    start = time.time()

    interpolator = FRIED_interpolator(logMstar=False, logF=False, logMdisk=False, logRdisk=False)

    end = time.time()

    print ('Building the interpolator for a {a} point grid took {b} s'.format(a=interpolator._grid.shape[0], b=end-start))


    start = time.time()

    Mdot = interpolator.interp(Mstar, F, Mdisk, Rdisk)

    end = time.time()

    print ('Performing {a} interpolations took {b} s, \nor {c} s per interpolation'.format(a=N, b=end-start, c=(end-start)/N))


    on_grid_mask = (interpolator._grid[:,0] == Mstar[0])*\
                   (interpolator._grid[:,1] == F[0])*\
                   (interpolator._grid[:,3] == Rdisk[0])

    Mdisk_on_grid = numpy.log10(interpolator._grid[ on_grid_mask, 2 ])


    plt.figure()
    plt.scatter(Mdisk_on_grid, interpolator._logMdot_grid[ on_grid_mask ])
    plt.plot(numpy.log10(Mdisk), numpy.log10(Mdot))

    plt.xlabel('log10 Disk Mass (MJup)')
    plt.ylabel('log10 Disk Mass Loss Rate (MSun/yr)')

    plt.title('Interpolations along disk mass slice, in lin space')

    plt.show()

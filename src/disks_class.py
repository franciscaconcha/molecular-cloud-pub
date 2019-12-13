from amuse.community.vader.interface import vader
from amuse.units import units, constants

import numpy


class Disk:

    def __init__(self,
                 disk_radius,
                 disk_mass,
                 alpha,
                 r_min=0.05 | units.au,
                 r_max=2000 | units.au,
                 n_cells=100,
                 linear=True):
        """ Initialize vader code for given parameters.

        :param disk_radius: disk radius. Must have units.au
        :param disk_mass: disk mass. Must have units.MSun
        :param alpha: turbulence parameter for viscosity, adimensional
        :param r_min: minimum radius of vader grid. Must have units.au
        :param r_max: maximum radius of vader grid. Must have units.au
        :param n_cells: number of cells for vader grid
        :param linear: linear interpolation
        :return: instance of vader code
        """
        self.disk = vader(redirection='none')
        self.disk.initialize_code()
        self.disk.initialize_keplerian_grid(n_cells,  # Number of cells
                                            linear,  # Linear?
                                            r_min,  # Grid Rmin
                                            r_max,  # Grid Rmax
                                            disk_mass  # Disk mass
                                            )

        sigma = self.column_density(disk_radius, disk_mass)
        self.disk.grid.column_density = sigma

        # The pressure follows the ideal gas law with a mean molecular weight of 2.33 hydrogen masses.
        # Since gamma ~ 1, this should not influence the simulation
        mH = 1.008 * constants.u
        T = 100. | units.K
        self.disk.grid.pressure = sigma * constants.kB * T / (2.33 * mH)

        self.disk.parameters.inner_pressure_boundary_type = 1
        self.disk.parameters.inner_pressure_boundary_torque = 0.0 | units.g * units.cm ** 2 / units.s ** 2
        self.disk.parameters.alpha = alpha
        self.disk.parameters.maximum_tolerated_change = 1E99

    def column_density(self,
                       rc,
                       mass,
                       lower_density=1E-12 | units.g / units.cm**2):
        """ Disk column density definition as in Eqs. 1, 2, and 3 of the paper.
            (Lynden-Bell & Pringle, 1974: Anderson et al. 2013)

        :param grid: disk grid
        :param rc: characteristic disk radius
        :param mass: disk mass
        :param lower_density: density limit for defining disk edge
        :return: disk column density in g / cm**2
        """
        r = self.disk.grid.r.value_in(units.au) | units.au
        rd = rc  # Anderson et al. 2013
        Md = mass

        Sigma_0 = Md / (2 * numpy.pi * rc ** 2 * (1 - numpy.exp(-rd / rc)))
        Sigma = Sigma_0 * (rc / r) * numpy.exp(-r / rc) * (r <= rc) + lower_density
        return Sigma

    def radius(self,
               density_limit=1E-10):
        """ Calculate the radius of a disk in a vader grid.

        :param disk: vader disk
        :param density_limit: density limit to designate disk border
        :return: disk radius in units.au
        """
        prev_r = self.disk.grid[0].r

        for i in range(len(self.disk.grid.r)):
            cell_density = self.disk.grid[i].column_density.value_in(units.g / units.cm ** 2)
            if cell_density < density_limit:
                return prev_r.value_in(units.au) | units.au
            prev_r = self.disk.grid[i].r

        return prev_r.value_in(units.au) | units.au

    def mass(self,
             radius):
        """ Calculate the mass of a vader disk inside a certain radius.

        :param disk: vader disk
        :param radius: disk radius to consider for mass calculation
        :return: disk mass in units.MJupiter
        """
        mass_cells = self.disk.grid.r[self.disk.grid.r <= radius]
        total_mass = 0

        for m, d, a in zip(mass_cells, self.disk.grid.column_density, self.disk.grid.area):
            total_mass += d.value_in(units.MJupiter / units.cm**2) * a.value_in(units.cm**2)

        return total_mass | units.MJupiter

    def density(self):
        """ Calculate the mean density of the disk, not considering the outer, low density limit.

        :param disk: vader disk
        :return: mean disk density in g / cm**2
        """
        radius = self.disk.radius()
        radius_index = numpy.where(self.disk.grid.r.value_in(units.au) == radius.value_in(units.au))
        density = self.disk.grid[:radius_index[0][0]].column_density.value_in(units.g / units.cm**2)
        return numpy.mean(density) | (units.g / units.cm**2)
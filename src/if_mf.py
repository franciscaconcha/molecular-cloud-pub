import numpy as np


def ionizing_flux (M, solar=True):
    '''
    Compute the ionizing flux of stars as a function of their mass, following Avedisova 1979
    Values are computed using linear interpolation in log space between data points

    M: masses of stars to compute ionizing flux of, in Solar masses (1D float array, shape (N))
    solar: use data for z=0.02 if True, use data for z=0.04 if False

    returns a 1D float array, shape (N), containing the ionising fluxes, in s^-1, corresponding to each stellar mass
    '''

    if solar:
        mass, logflux = np.loadtxt('ionising_flux_z02.txt', unpack=True)
    else:
        mass, logflux = np.loadtxt('ionising_flux_z04.txt', unpack=True)

    mass = mass[::-1]
    logflux = logflux[::-1]

    logmass = np.log10(mass)

    N = len(M)
    n = len(mass)

    Ndot = np.zeros(N)

    for i in range(n-1):

        mask = (M >= mass[i])*(M < mass[i+1])

        #print (np.sum(mask))

        Ndot[ mask ] = 10.**( (logflux[i+1] - logflux[i])/(logmass[i+1] - logmass[i])*(np.log10(M[ mask ]) - logmass[i]) + logflux[i] )

    return Ndot


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    M = np.logspace(0., 2., num=1000)

    Ndot = ionizing_flux(M)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    ax.plot(M[ Ndot > 0. ], Ndot[ Ndot > 0. ])

    mass, logflux = np.loadtxt('ionising_flux_z02.txt', unpack=True)

    ax.scatter(mass, 10.**logflux)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Stellar Mass [MSun]')
    ax.set_ylabel('Ionising Flux [s$^{-1}$]')

    plt.show()

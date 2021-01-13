# Evolution of circumstellar discs in young star forming regions #
#### Francisca Concha-Ramírez, Simon Portegies Zwart, Martijn J. C. Wilhelm
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3897171.svg)](https://doi.org/10.5281/zenodo.3897171)
 [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) ![python](https://img.shields.io/badge/python-3.0-yellow.svg)


Project script for the paper: Evolution of circumstellar discs in young star forming regions (link TBD)

These scripts simulate the collapse of a giant molecular cloud with star formation, and evolves circumstellar discs around the stars. The discs are subject to viscous expansion, dynamical encounters, external and internal photoevaporation, and internal dust evolution.

The simulations are implemented in three scripts to be run sequentially. Below we explain how to set up and run each of the scripts.

For questions please contact Francisca Concha-Ramírez, e-mail: hello at francisca.cr

## Getting Started

### Requisites
* Python 2.7. Should work fine with Python 3 but it has not been tested.
* AMUSE: https://amusecode.github.io
* VADER: https://bitbucket.org/krumholz/vader/src
* scipy

### How to set up

You need to install AMUSE with the developer option so that you can access the source code.
Then, download VADER and put it inside the folder: ```/amuse/src/amuse/community```. 
In the ```/vader``` folder in this repository there are several files related to VADER. They should go in the following directories:

* ```interface.cc``` and ```interface.py``` should go on ```/amuse/src/amuse/community/vader/```
* ```userFunc_pedisk.c``` and ```userFunc_none.c``` should go on ```amuse/src/amuse/community/vader/src/prob``` (these 2 files are redundant but it was a way to go around some compilation issues)
* ```Makefile_interface``` should go on ```/amuse/src/amuse/community/vader/``` and renamed ```Makefile```
* ```Makefile_source``` should go on ```/amuse/src/amuse/community/vader/src/``` and renamed ```Makefile```
* Compile VADER from the main AMUSE directory with ```make vader.code```
* Add this line to ```amuse/src/amuse/lab.py```: 
```
from amuse.community.vader.interface import vader, vaderInterface
```

You should now be ready to run the simulation scripts.

### Running the simulations

#### 1. Molecular cloud collapse

The first step of the simulations corresponds to the collapse of a giant molecular cloud. This step is implemented in `src/molecular_cloud_collapse_with_stars.py`. You can run an individual simulation by using the AMUSE script directly from the home directory:

```
amuse.sh src/molecular_cloud_collapse_with_stars.py
```

The script has extensive options which can be passed through the command line. For a full list of these options run: `amuse.sh src/molecular_cloud_collapse_with_stars.py --help`.

* `-f`: filename to use to restart the runs
* `-s`: path to save the results
* `--tend`: end time for the simulation (default 5 Myr)
* `--dt_diag`: dt to print diagnostics (default 0.1 Myr)
* `--Ncloud`: number of SPH particles (default 1000)
* `--Mcloud`: mass of the cloud, in solar masses (default 1E3 MSun)
* `--Rcloud`: initial radius of the cloud, in parsec (default 3 parsec)
* `--fcloud`: factor within which to locate the newly formed stars. When a star forms from a sink, its initial location will be randomly determined within a sphere of radius `fcloud * sink_radius` (default 1)

This script will output three types of files: `hydro_gas_particles_i*.amuse` (gas particles), `hydro_sink_particles_i*.amuse` (sink particles), and `hydro_stars_particles_i*.amuse` (stars) where `i` is followed by an index. The script stops when all the stars in the IMF have been formed (see paper for details) or when `tend` is reached, whichever happens first. 

#### 2. Evolving the discs during star formation

The first script deals only with the collapse of the molecular cloud and star formation, and with the dynamics of the newly formed stars. To evolve the circumstellar discs during this time frame, we implemented a separate script in `src/presteps.py`. This script reads the `hydro_stars_particles_i*.amuse` files created by the previous step, reads stars as they are born, assigns them circumstellar discs and evolves them in time. This script can be run as:

```amuse.sh src/presteps.py```

The script has extensive options which can be passed through the command line. For a full list of these options run: `amuse.sh src/presteps.py --help`. The main options are:

* `-p`: path to the `hydro_stars_particles_i*.amuse` files to read
* `-s`: path to save the results
* `-f`: path to the FRIED grid files, for external photoevaporation ([source](https://arxiv.org/abs/1808.07484)) 
* `-c`: number of cores to use (default 2)  

This script will output a series of files named `N<final number of stars>_t<timestamp>Myr.hdf5`. Each of these files contains the dynamics and disc information for all stars formed before the timestamp. This script ends at the same timestamp that the molecular cloud collapse script ends.

#### 3. Evolution after star formation has ended

Once star formation ends, the gas in the cloud is removed instantaneously and the stars and discs continue evolving for 2 Myr. This third step of the simulations is implemented in the script `src/vader_cluster_parallel.py`. This script will take as input the last `N<final number of stars>_t<timestamp>Myr.hdf5` file generated in the last step and will continue the evolution of the stellar dynamics, stellar evolution, disc expansion, internal and external photoevaporation, and dust inside the disc.

This script takes a series of options which can be listed with `amuse.sh src/vader_cluster_parallel.py --help`. The main options are:

* `-p`: path to the `N<final number of stars>_t<timestamp>Myr.hdf5` files to read
* `-s`: path to save the results
* `-f`: path to the FRIED grid files, for external photoevaporation ([source](https://arxiv.org/abs/1808.07484)) 
* `-c`: number of cores to use (default 2)  


## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details

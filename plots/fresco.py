import matplotlib.pyplot as plt
from amuse.datamodel import Particles
from amuse.units import units, nbody_system
from amuse.community.sse.interface import SSE
from amuse.ext.masc import make_a_star_cluster
from amuse.ext.fresco import make_fresco_image
#?make_fresco_image  # See options
stars_file = 'results/07082020/0/gravity_stars_t0.01Myr.hdf5'
stars = read_set_from_file(stars_file, 'hdf5', close_file=True)
gas = Particles()
se = SSE()
se.particles.add_particles(stars)
from_se = se.particles.new_channel_to(stars)
from_se.copy()
image, vmax = make_fresco_image(
    stars, gas,
    mode=["stars"],
    return_vmax=True,
)
plt.imshow(image)
plt.show()
from amuse.lab import *

prev = 0 | units.Myr

#for i in range(13, 31):
i=20
stars = read_set_from_file("hydrostars/1/hydro_stars_particles_i00{0}.amuse".format(i), "hdf5", close_file=True)

# Set particle masses
stars[stars.stellar_mass <= 1.9 | units.MSun].disk_mass = 0.1 * stars[stars.stellar_mass <= 1.9 | units.MSun].stellar_mass
stars[stars.stellar_mass > 1.9 | units.MSun].disk_mass = 0.0 | units.MSun
stars.mass = stars.stellar_mass + stars.disk_mass

# Set disk radii
stars[stars.stellar_mass <= 1.9 | units.MSun].disk_radius = 100 | units.au
stars[stars.stellar_mass > 1.9 | units.MSun].disk_radius = 0.0 | units.au


t_end = 2 | units.Myr
t = 0 | units.Myr
dt = 1000 | units.yr

tmin = 20 | units.Myr
for s in stars:
    if s.tborn < tmin:
        tmin = s.tborn

print("first stars at ", tmin.in_(units.Myr))
first_stars = stars[stars.tborn == tmin]

converter = nbody_system.nbody_to_si(stars.stellar_mass.sum(), 0.5 | units.parsec)

gravity = ph4(converter)
gravity.parameters.timestep_parameter = 0.01
gravity.parameters.epsilon_squared = (100 | units.au) ** 2
gravity.particles.add_particles(first_stars)

channel_from_gravity_to_framework = gravity.particles.new_channel_to(stars)
channel_from_framework_to_gravity = stars.new_channel_to(gravity.particles,
                                                         attributes=['collisional_radius'],
                                                         target_names=['radius'])

tprev = tmin  # Otherwise I will add first_stars twice

while t < t_end:
    print("t = {0}".format(t.in_(units.Myr)))
    if t > tmin:
        prev_stars = stars[stars.tborn > tprev]
        new_stars = prev_stars[prev_stars.tborn <= t]
        gravity.particles.add_particles(new_stars)
        channel_from_gravity_to_framework.copy()
        tprev = t

    gravity.evolve_model(t + dt)

    print(t, len(gravity.particles))
    t += dt

from amuse.lab import *
import matplotlib.pyplot as pyplot
import numpy

prev = 0 | units.Myr

colors = ['#1373b2', '#42a6cc', '#7cccc4', '#b4e2ba', '#daf0d4', 'k']

decay = []
times = []

prev_keys = None

for i in range(14, 28):
    sinks = read_set_from_file("results/23072020/0/hydro_sink_particles_i00{0}.amuse".format(i),
                               "hdf5", close_file=True)

    sink = sinks[sinks.key == 13317617281889275833]
    #print sink

    """if i == 14:
        prev_keys = sinks.key

    print i
    print list(set(prev_keys).intersection(sinks.key))"""

    time = sinks.get_timestamp().value_in(units.Myr)
    #sink = sinks[sinks.key == sink_key]

    if sink.time_threshold:
        decay.append(sink.time_threshold.value_in(units.Myr) - time)
        times.append(sinks.get_timestamp().value_in(units.Myr))
    else:
        decay.append(0.0)
        times.append(sinks.get_timestamp().value_in(units.Myr))

    print time, sink.merged_keys

print decay
print len(decay), len(times[:len(decay)])

pyplot.plot(times[:len(decay)], decay)
pyplot.show()
import numpy
import matplotlib.pyplot as plt
import sys

sys.path.insert(0,'..')

import vtr_to_numpy as vtn

def mean_velx_time(name):
    fig = plt.figure(figsize=(8,8))
    data = vtn.vtrToNumpy(vtn.loadFromResults(), name)[0]
    arr = numpy.average(data, axis = (0,1))
    fig.plot(arr)
    fig.show()

mean_velx_time("Electron - Solar wind-temperature")

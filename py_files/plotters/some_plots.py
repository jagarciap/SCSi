import numpy
import matplotlib.pyplot as plt
import sys

sys.path.insert(0,'..')

import constants as c
from mesh import Mesh_2D_rm
import vtr_to_numpy as vtn

mesh = Mesh_2D_rm(c.XMIN, c.XMAX, c.YMIN, c.YMAX, c.NX, c.NY, c.DEPTH)

def mean_velx_time(name):
    fig = plt.figure(figsize=(8,8))
    data = vtn.vtrToNumpy(mesh, vtn.loadFromResults(), name)
    arr_e = numpy.average(data[0], axis = 0)
    arr_p = numpy.average(data[1], axis = 0)
    plt.plot(arr_e, label = 'Electrons')
    plt.plot(arr_p, label = 'Protons')
    plt.legend()
    fig.show()

mean_velx_time(["Electron - Solar wind-temperature", "Proton - Solar wind-temperature"])

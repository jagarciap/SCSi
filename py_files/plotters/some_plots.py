import numpy
import matplotlib.pyplot as plt
import pdb
import sys

sys.path.insert(0,'..')

import constants as c
from mesh import Mesh_2D_rm
import vtr_to_numpy as vtn

mesh = Mesh_2D_rm(c.XMIN, c.XMAX, c.YMIN, c.YMAX, c.NX, c.NY, c.DEPTH)

def mean_temperature_time(name):
    fig = plt.figure(figsize=(8,8))
    data = vtn.vtrToNumpy(mesh, vtn.loadFromResults(), name)
    arr_e = numpy.average(data[0], axis = 0)
    arr_p = numpy.average(data[1], axis = 0)
    plt.plot(arr_e, label = 'Electrons')
    plt.plot(arr_p, label = 'Protons')
    plt.legend()
    plt.show()

def mean_vel_time(name):
    fig = plt.figure(figsize=(8,8))
    data = vtn.vtrToNumpy(mesh, vtn.loadFromResults(), name)
    arr_e = numpy.average(numpy.linalg.norm(data[0], axis = 1), axis = 0)
    arr_p = numpy.average(numpy.linalg.norm(data[1], axis = 1), axis = 0)
    plt.plot(arr_e, label = 'Electrons')
    plt.plot(arr_p, label = 'Protons')
    plt.legend()
    plt.show()

def mean_density_time(name):
    fig = plt.figure(figsize=(8,8))
    data = vtn.vtrToNumpy(mesh, vtn.loadFromResults(), name)
    arr_e = numpy.average(data[0], axis = 0)
    arr_p = numpy.average(data[1], axis = 0)
    plt.plot(arr_e, label = 'Electrons')
    plt.plot(arr_p, label = 'Protons')
    plt.legend()
    plt.show()

def extremes_rho_time(name):
    fig = plt.figure(figsize=(8,8))
    data = vtn.vtrToNumpy(mesh, vtn.loadFromResults(), name)
    arr_min = numpy.min(data[1]-data[0], axis = 0)
    arr_max = numpy.max(data[1]-data[0], axis = 0)
    plt.plot(arr_min, label = 'min')
    plt.plot(arr_max, label = 'max')
    plt.legend()
    plt.show()

def debye_length_test(name):
    fig = plt.figure(figsize=(8,8))
    data = vtn.vtrToNumpy(mesh, vtn.loadFromResults(), name)
    left = c.NX*(int(c.NY/2)+1)
    cut = data[0][left:left+c.NX,-1]
    x = numpy.linspace(c.XMIN, c.XMAX, num = c.NX)
    def debye_function(x):
        lambda_d = 1/numpy.sqrt(c.E_N*c.QE*c.QE/c.EPS_0/c.K/c.E_T+c.P_N*c.QE*c.QE/c.EPS_0/c.K/c.P_T)
        q0 = -c.QE*1.6e12*(c.YMAX-c.YMIN)/(c.NY-1)*(c.XMAX-c.XMIN)/(c.NX-1)*c.DEPTH
        return 1/4/numpy.pi/c.EPS_0*q0/abs(x)*numpy.exp(-abs(x)/lambda_d)*1/2.1
    theory = debye_function(x-(c.XMAX-c.XMIN)/2)
    plt.plot(x, theory, color = 'black', label = 'Theory')
    plt.scatter(x, cut, color = 'red', marker = '.', label = 'Simulation')
    plt.legend()
    plt.show()


#extremes_rho_time(["Electron - Solar wind-density", "Proton - Solar wind-density"])
#mean_density_time(["Electron - Solar wind-density", "Proton - Solar wind-density"])
#mean_temperature_time(["Electron - Solar wind-temperature", "Proton - Solar wind-temperature"])
#mean_vel_time(["Electron - Solar wind-velocity", "Proton - Solar wind-velocity"])
debye_length_test(["Electric - Electrostatic_2D_rm-potential"])

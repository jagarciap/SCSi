import numpy
import matplotlib.pyplot as plt
import pdb
from scipy import stats, optimize
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

def debye_function(x, *args):
    q0 = -c.QE*1.6e12*(c.YMAX-c.YMIN)/(c.NY-1)*(c.XMAX-c.XMIN)/(c.NX-1)*c.DEPTH
    return args[0]*1/4/numpy.pi/c.EPS_0*q0/abs(x)*numpy.exp(-abs(x)/args[1])

def debye_shielding_fit(x, y, guess):
    return optimize.curve_fit(debye_function,x,y, p0=guess)

def debye_length_test(name):
    fig = plt.figure(figsize=(8,8))
    #Preparing data
    data = vtn.vtrToNumpy(mesh, vtn.loadFromResults(), name)
    left = c.NX*(int(c.NY/2)+1)
    cut = data[0][left:left+c.NX,-1]
    x = numpy.linspace(c.XMIN, c.XMAX, num = c.NX)
    offset = (c.XMAX-c.XMIN)/2
    dx = (c.XMAX-c.XMIN)/(c.NX-1)
    #Omitting the center
    ind = numpy.where(numpy.logical_or(x < offset-4*dx, x > offset+4*dx))

    #Debye Lenght study
    lambda_d = 1/numpy.sqrt(c.E_N*c.QE*c.QE/c.EPS_0/c.K/c.E_T+c.P_N*c.QE*c.QE/c.EPS_0/c.K/c.P_T)
    print("Theory", lambda_d)
    theory = debye_function(x-offset, 1.0, lambda_d)
    guess = (1.0, lambda_d)
    params, errors = debye_shielding_fit(x[ind]-offset, cut[ind], guess)

    plt.plot(x, theory, color = 'black', label = 'Theory')
    plt.scatter(x, cut, color = 'red', marker = '.', label = 'Simulation')
    plt.plot(x, debye_function(x-offset, *params), color = 'blue', label = 'Fit')
    plt.legend()
    plt.show()

    #Further analysis
    fig = plt.figure(figsize=(8,8))
    val = []
    err = []
    for i in range(1,10):
        ind = numpy.where(numpy.logical_or(x < offset-i*dx, x > offset+i*dx))
        guess = (1.0, lambda_d)
        params, errors = debye_shielding_fit(x[ind]-offset, cut[ind], guess)
        val.append(params[1])
        err.append(numpy.sqrt(numpy.diag(errors)[1]))
        print(params, numpy.sqrt(numpy.diag(errors)))
    plt.axhline(y=lambda_d, color = 'black')
    plt.errorbar(numpy.arange(1,10), val, yerr=err, marker = '.', linestyle = '', ecolor = 'black')
    plt.show()


extremes_rho_time(["Electron - Solar wind-density", "Proton - Solar wind-density"])
mean_density_time(["Electron - Solar wind-density", "Proton - Solar wind-density"])
mean_temperature_time(["Electron - Solar wind-temperature", "Proton - Solar wind-temperature"])
mean_vel_time(["Electron - Solar wind-velocity", "Proton - Solar wind-velocity"])
#debye_length_test(["Electric - Electrostatic_2D_rm-potential"])

#This file contains the functions used to plot the files in "results_particles"
import numpy
import matplotlib.pyplot as plt
import os
import pdb
from scipy import stats, optimize
from subprocess import check_output
import sys

import constants as c

## ---------------------------------------------------------------------------------------------------------------
# Functions for handling Gaussian curves.
# SOURCE: 'DataAnalysis2.py' in '/home/jorge/Documents/Instituto_Balseiro/Semestre2017-1/Experimental_IV/Mediciones_en_coincidencia'.
#           There are addiional functions there that can be of use in the future.
## ---------------------------------------------------------------------------------------------------------------

def nGaussianFunc(x, *params):
    y = numpy.zeros_like(x,dtype=numpy.float)
    for i in range (0, len(params), 3):
        #print(params[i], params[i+1], params[i+2])
        y += params[i]*numpy.exp(-(x-params[i+1])**2/(2*params[i+2]**2))
    return y

def nGaussianAreas(*params):
    y = []
    for i in range(0,len(params),3):
        y.append(params[i]*numpy.sqrt(2*numpy.pi)*params[i+2])
    return y

def nGaussianFit (x, y, guess):
    return optimize.curve_fit(nGaussianFunc,x,y, p0=guess)

## ---------------------------------------------------------------------------------------------------------------

def loadFromResults():
    cwd_base = os.getcwd().rsplit(sep = os.sep, maxsplit = 1)
    cwd = os.path.join(cwd_base[0], 'results_particles','')
    stdout = check_output('ls' +' {}'.format(cwd), shell=True)
    files = stdout.decode().split(sep='\n')
    for i in range(len(files)):
        files[i] = os.path.join(cwd, files[i])
    return files[:-1]

def phase_space(filename):
    array = numpy.loadtxt(filename)
    fig = plt.figure(figsize=(8,5))
    #plt.scatter(array[:,0]+array[:,1], array[:,3], marker = '.', label = 'electrons') 
    plt.scatter(array[:,4]+array[:,5], array[:,7], marker = '.', label = 'protons') 
    plt.legend()
    plt.title('y_Vel-Pos phase space')
    plt.show()
    fig.savefig('phase_space_vely_pos_2.png')

def vel_distribution_electrons(filename):
    array = numpy.loadtxt(filename)
    fig = plt.figure(figsize=(8,5))
    #datax = plt.hist(array[:,2],81, alpha = 0.5, label = 'electrons, vel_x') 
    #datay = plt.hist(array[:,3],81, alpha = 0.5, label = 'electrons, vel_y') 
    datamag = plt.hist(numpy.sqrt(array[:,2]*array[:,2]+array[:,3]*array[:,3]), 41, alpha = 0.5, label = 'electrons, vel_mag')
    ##Gaussian fits
    #params, errors = nGaussianFit(datax[1][:-1], datax[0], [datax[0].max(), 3e5, numpy.sqrt(c.K*c.E_T/c.ME)])
    #x = numpy.linspace(datax[1].min(), datax[1].max(), num = 100)
    #plt.plot(x, nGaussianFunc(x, *params,), label = 'velx_fit')
    #print("electrons_vel_x")
    #print(params, errors)
    #print("Temperature_vel_x", params[2]**2*c.ME/c.K/c.EV_TO_K)
    #params, errors = nGaussianFit(datay[1][:-1], datay[0], [datay[0].max(), 0, numpy.sqrt(c.K*c.E_T/c.ME)])
    #x = numpy.linspace(datay[1].min(), datay[1].max(), num = 100)
    #plt.plot(x, nGaussianFunc(x, *params,), label = 'vely_fit')
    #print("electrons_vel_y")
    #print(params, errors)
    #print("Temperature_vel_y", params[2]**2*c.ME/c.K/c.EV_TO_K)

    plt.legend()
    plt.title('vel_distribution_electrons')
    plt.show()
    #fig.savefig('vel_distribution_electrons.png')

def vel_distribution_protons(filename):
    array = numpy.loadtxt(filename)
    fig = plt.figure(figsize=(8,5))
    #datax = plt.hist(array[:,6],81, alpha = 0.5, label = 'protons, vel_x') 
    #datay = plt.hist(array[:,7],81, alpha = 0.5, label = 'protons, vel_y') 
    datamag = plt.hist(numpy.sqrt(array[:,6]*array[:,6]+array[:,7]*array[:,7]), 81, alpha = 0.5, label = 'protons, vel_mag')
    ##Gaussian fits
    #params, errors = nGaussianFit(datax[1][:-1], datax[0], [datax[0].max(), 3e5, numpy.sqrt(c.K*c.P_T/c.MP)])
    #x = numpy.linspace(datax[1].min(), datax[1].max(), num = 100)
    #plt.plot(x, nGaussianFunc(x, *params,), label = 'vel_x_fit')
    #print("protons_vel_x")
    #print(params, errors)
    #print("Temperature_vel_x", params[2]**2*c.MP/c.K/c.EV_TO_K)
    #params, errors = nGaussianFit(datay[1][:-1], datay[0], [datay[0].max(), 0, numpy.sqrt(c.K*c.P_T/c.MP)])
    #x = numpy.linspace(datay[1].min(), datay[1].max(), num = 100)
    #plt.plot(x, nGaussianFunc(x, *params,), label = 'vel_y_fit')
    #print("protons_vel_y")
    #print(params, errors)
    #print("Temperature_vel_y", params[2]**2*c.MP/c.K/c.EV_TO_K)

    plt.legend()
    plt.title('vel_distribution_protons')
    plt.show()
    #fig.savefig('vel_distribution_protons.png')

vel_distribution_electrons(loadFromResults()[-1])
vel_distribution_protons(loadFromResults()[-1])

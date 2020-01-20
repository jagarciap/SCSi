# File in charge of printing
import os
import copy
from datetime import datetime
import motion
import pyevtk.hl as vtk
import numpy
import pickle
import constants as c
import pdb


#NOTE: to keep the good detachment and modularity of the program, this function can be generalized easily. Each class that has something to print has a print function, which returns a dictionary.
#       then, here, I just need to unravel and packed everything again.
def output_vtk(ts, mesh, electrons, protons, e_field):
    #Creating domain
    nx = mesh.nx
    ny = mesh.ny
    i = numpy.arange(0, nx, dtype ='int16')
    j = numpy.arange(0, ny, dtype ='int16')
    temp = numpy.zeros((1), dtype = 'int16')

    #Files created
    cwd = os.path.split(os.getcwd())[0]
    vtkstring = cwd+'/results/ts{:d}'.format(ts)
    vtk.gridToVTK(vtkstring, i, j, temp, pointData = {
        'n_e': electrons.mesh_values.density.reshape((nx, ny, 1), order = 'F'),\
        'e_vel': (numpy.reshape(copy.copy(electrons.mesh_values.velocity[:,0]),(nx,ny,1), order = 'F'),\
                numpy.reshape(copy.copy(electrons.mesh_values.velocity[:,1]),(nx,ny,1), order = 'F'), \
                numpy.zeros((nx,ny,1))),\
        'n_p': protons.mesh_values.density.reshape((nx, ny, 1), order = 'F'),\
        'p_vel': (numpy.reshape(copy.copy(protons.mesh_values.velocity[:,0]),(nx,ny,1), order = 'F'),\
                numpy.reshape(copy.copy(protons.mesh_values.velocity[:,1]),(nx,ny,1), order = 'F'), \
                numpy.zeros((nx,ny,1))),\
        'e_field': (numpy.reshape(copy.copy(e_field.field[:,0]),(nx,ny,1), order = 'F'),\
                numpy.reshape(copy.copy(e_field.field[:,1]),(nx,ny,1), order = 'F'), \
                numpy.zeros((nx,ny,1)))})

# The function prints a file for a particular timestep 'ts' where the species being tracked are printed. Columns are for each component of each species, so for 2D:
#   specie1.x \t specie1.y \t specie2.x etc. Each row is a different particle for a particular species.
def particle_tracker(ts, *args):
    # Checking tracking method
    for spc in args:
        if spc.part_values.current_n > spc.part_values.num_tracked and numpy.any(spc.part_values.trackers == spc.part_values.max_n):
            print("Error in species: ", spc.type)
            print(spc.part_values.current_n, spc.part_values.num_tracked)
            pdb.set_trace()
            raise ValueError("There should not be any invalid values")

    #Creating array to be printed and the header
    narray = numpy.zeros((args[0].part_values.num_tracked, args[0].pos_dim*len(args)))
    nHeader = ''
    for i in range(len(args)):
        ind = numpy.argwhere(args[i].part_values.trackers != args[i].part_values.max_n)
        narray[ind, args[i].pos_dim*i:args[i].pos_dim*(i+1)] = args[i].part_values.position[args[i].part_values.trackers[ind],:]
        nHeader += args[i].type + '\t'

    cwd = os.path.split(os.getcwd())[0]
    workfile = cwd+'/particle_tracker/ts={:05d}.dat'.format(ts)
    nHeader = 'No. of particles = {:d} \n'.format(args[0].part_values.num_tracked)+nHeader
    numpy.savetxt(workfile, narray , fmt = '%.5e', delimiter = '\t', header = nHeader)

def save_current_state(ts, world, ions, electrons):
    cwd = os.getcwd()
    time = datetime.now().strftime('%Y-%m-%d_%Hh%Mm')
    string = cwd+'/previous_executions/sys_ts={:d}_'.format(ts)+time+'.pkl'
    with open (string, 'wb') as output:
        pickle.dump(ts, output, -1)
        pickle.dump(world, output, -1)
        pickle.dump(ions, output, -1)
        pickle.dump(electrons, output, -1)

def load_state(filename):
    cwd = os.getcwd()
    filename = cwd+'/previous_executions/'+filename
    with open (filename, 'rb') as pinput:
        ts = pickle.load(pinput)
        world = pickle.load(pinput)
        ions = pickle.load(pinput)
        electrons = pickle.load(pinput)
        return ts, world, ions, electrons


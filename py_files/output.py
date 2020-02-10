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
def saveVTK(ts, mesh, electrons, protons, e_field):
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

#       +Method that loads the information of the system from a '.vtk' and stores it in the arguments *args.
#       +Structure to be followed in *args:
#       ++ts (timestep); fields: Electrics, Magnetics; Species: Electrons, Protons, Ions, Neutrals.
#       ++Inside the types not further specified now, an alphabetical order with respect to the classes' names will be maintained.
def loadVTK(filename, *args):
    #Preparing path
    cwd = os.path.split(os.getcwd())[0]
    filename = cwd+'/initial_conditions/'+filename


# The function prints a file for a particular timestep 'ts' where the species being tracked are printed. Columns are for each component of each species, so for 2D:
#   specie1.x \t specie1.y \t specie2.x etc. Each row is a different particle for a particular species.
def particleTracker(ts, *args):
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

#       +Method that stores the information of the system, given in *args, in a '.pkl' file. See 'Pickle' module for further information.
#       +Structure to be followed in *args:
#       ++ts (timestep); fields: Electrics, Magnetics; Species: Electrons, Protons, Ions, Neutrals; part_solver (Particle Solver).
#       ++Inside the types not further specified now, an alphabetical order with respect to the classes' names will be maintained.
def savePickle(*args):
    #Creating file's name
    cwd = os.path.split(os.getcwd())[0]
    time = datetime.now().strftime('%Y-%m-%d_%Hh%Mm')
    string = cwd+'/previous_executions/sys_ts={:d}_'.format(ts)+time+'.pkl'
    #Storing information
    with open (string, 'wb') as output:
        for arg in args:
            pickle.dump(arg, output, -1)

#       +Method that loads the information of the system from a '.pkl' and stores it in the arguments *args. See 'Pickle' module for further information.
#       +Structure to be followed in *args:
#       ++ts (timestep); fields: Electrics, Magnetics; Species: Electrons, Protons, Ions, Neutrals; part_solver (Particle Solver).
#       ++Inside the types not further specified now, an alphabetical order with respect to the classes' names will be maintained.
def loadPickle(filename, *args):
    #Preparing path
    cwd = os.path.split(os.getcwd())[0]
    filename = cwd+'/initial_conditions/'+filename
    with open (filename, 'rb') as pinput:
        for arg in args:
            arg = pickle.load(pinput)

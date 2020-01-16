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

def particle_Tracker_init(ions, electrons, num = 100):
    ions_step = ions.np//num
    electrons_step = electrons.np//num
    ions.part.printable[numpy.arange(0,ions_step*num,ions_step)] = numpy.arange(1,num+1)
    electrons.part.printable[numpy.arange(0, electrons_step*num, electrons_step)] = numpy.arange(1,num+1)

def particle_Tracker_print(ts, ions, neutrals, num = 100):
    
    ind_ions = numpy.flatnonzero(ions.part.printable[:ions.np])
    rem_ions = num-len(ind_ions)
    ind_neutrals = numpy.flatnonzero(neutrals.part.printable[:neutrals.np])
    rem_neutrals = num-len(ind_neutrals)

    a, test_ions = numpy.unique(ions.part.print_ind, return_counts = True)
    bool_ions = numpy.any(test_ions > 1)
    a, test_neutrals = numpy.unique(neutrals.part.print_ind, return_counts = True)
    bool_neutrals = numpy.any(test_neutrals > 1)
    if rem_neutrals != len(neutrals.part.print_ind) or rem_ions != len(ions.part.print_ind) \
            or bool_ions or bool_neutrals:
                pdb.set_trace()

    if rem_ions > 0:
        step = (ions.np-ind_ions[-1])//rem_ions
        for i in range(rem_ions):
            ions.part.printable[ind_ions[-1]+1+i*step] = ions.part.print_ind.pop(0)

    if rem_neutrals > 0:
        step = (neutrals.np-ind_neutrals[-1])//rem_neutrals
        for i in range(rem_neutrals):
            neutrals.part.printable[ind_neutrals[-1]+1+i*step] = neutrals.part.print_ind.pop(0)

    ind_ions = numpy.flatnonzero(ions.part.printable[:ions.np])
    ind_neutrals = numpy.flatnonzero(neutrals.part.printable[:neutrals.np])
    np_ions = len(ind_ions)
    np_neutrals = len(ind_neutrals)
    if len(ions.part.print_ind) > 0 or len(neutrals.part.print_ind) > 0:
        pdb.set_trace()

    narray = numpy.zeros((num,4))
    narray[ions.part.printable[ind_ions]-1, 0] = ions.part.x[ind_ions,0]
    narray[ions.part.printable[ind_ions]-1, 1] = ions.part.x[ind_ions,1]
    narray[neutrals.part.printable[ind_neutrals]-1, 2] = neutrals.part.x[ind_neutrals,0]
    narray[neutrals.part.printable[ind_neutrals]-1, 3] = neutrals.part.x[ind_neutrals,1]

    cwd = os.getcwd()
    workfile = cwd+'/particle_tracker/ts={:05d}.dat'.format(ts)
    nheader = 'Max. number of particles: \t {:4d} \n Ions \t Neutrals \n {:3d} \t {:3d} \n'.format( num,np_ions, np_neutrals)
    numpy.savetxt(workfile, narray , fmt = '%.5e', delimiter = '\t', header = nheader)


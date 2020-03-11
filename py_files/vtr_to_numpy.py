# This file provides methods to transform .vtr files into numpy arrays
import numpy
import os
import pdb
from subprocess import check_output
import sys
from vtk.util.numpy_support import vtk_to_numpy

def vtrToNumpy(mesh, filenames, names):
    arrays = []
    #First filename
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    cwd = os.path.split(os.getcwd())[0]
    filename = cwd+'/results/'+filenames[0]
    reader = mesh.vtkReader()
    reader.SetFileName(filename)
    reader.Update()
    output = reader.GetOutput()
    for name in names:
        temp = mesh.reverseVTKOrdering(vtk_to_numpy(output.GetPointData().GetArray(name)))
        temp = numpy.expand_dims(temp, axis = temp.ndim)
        arrays.append(temp)
    for filename in filenames[1:]:
        cwd = os.path.split(os.getcwd())[0]
        filename = cwd+'/results/'+filename
        reader = mesh.vtkReader()
        reader.SetFileName(filename)
        reader.Update()
        output = reader.GetOutput()
        for i in range(len(names)):
            temp = mesh.reverseVTKOrdering(vtk_to_numpy(output.GetPointData().GetArray(names[i])))
            temp = numpy.expand_dims(temp, axis = temp.ndim)
            arrays[i] = numpy.append(arrays[i], temp, axis = arrays[i].ndim-1)
    return arrays

def loadFromResults():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    cwd = os.path.split(os.getcwd())[0]
    cwd = cwd+'/results/'
    stdout = check_output('ls' +' {}'.format(cwd), shell=True)
    files = stdout.decode().split(sep='\n')
    return files[:-1]


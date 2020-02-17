# This file provides methods to transform .vtr files into numpy arrays
import numpy
import os
import pdb
import sys
from vtk.util.numpy_support import vtk_to_numpy
from subprocess import check_output

from mesh import Mesh_2D_rm as mesh

def vtrToNumpy(filenames, names):
    arrays = []
    #First filename
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    cwd = os.path.split(os.getcwd())[0]
    filename = cwd+'/results/'+filenames[0]
    reader = mesh.vtkReader(mesh)
    reader.SetFileName(filename)
    reader.Update()
    output = reader.GetOutput()
    for name in names:
        arrays.append(mesh.reverseVTKOrdering(mesh,vtk_to_numpy(output.GetPointData().GetArray(self.name)))[:,None])
    for filename in filenames[1:]:
        cwd = os.path.split(os.getcwd())[0]
        filename = cwd+'/results/'+filename
        reader = mesh.vtkReader(mesh)
        reader.SetFileName(filename)
        reader.Update()
        output = reader.GetOutput()
        for i in range(len(names)):
            arrays[i] = numpy.append(arrays[i],mesh.reverseVTKOrdering(mesh,vtk_to_numpy(output.GetPointData().GetArray(self.names[i]))), axis = 2)
    return arrays

def loadFromResults():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    cwd = os.path.split(os.getcwd())[0]
    cwd = cwd+'/results/'
    stdout = check_output('ls' +' {}'.format(cwd), shell=True)
    files = stdout.decode().split(sep='\n')
    return files


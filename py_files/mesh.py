#Data structures to hold domain information
import copy
import numpy
import os
import pdb
import vtk

import evtk.hl as vtk

import Boundaries.outer_2D_rectangular as ob

#Mesh (Abstract)(Association between Mesh and PIC):
#
#Definition = Defines the type of mesh.
#Attributes:
#	+nPoints (int) = Number of points in the mesh.
#       +boundaries ([Boundary]) = List of different boundaries that define the mesh.
#	+volumes ([double]) = Volume of each node.
#Methods:
#       +setDomain() = This function, with the values provided by the boundary files, will create the mesh, by setting up volumes, nPoints, boundaries and any other subclass variable.
#	+getPosition([int] i): [double, double y] = For a each index return its real position.
#	+getIndex([double,double] pos): [double,double] = For each real position returns its index value. Important to remember that the number of columns may vary
#           depending on the actual type of mesh subclass used.
#	+arrayToIndex([ind] array): [int, int] = For the indexes in the 1D array, obtain the indexes used for the particular mesh.
#	+indexToArray([ind, ind] index): [int] = For the indexes used for the particular mesh, obtain the 1D version for the array.
#       +vtkOrdering(array): array = The array received as argument is ordered in such a way it can be stored ina VTK file.
#           The result is returned as a new array.
#       +vtkOrdering(array): array = The array received as argument comes with vtk ordering and is reshaped to be stored properly in the code.
#       +vtkReader(): Reader = Return the reader from module vtk that is able to read the vtk file.
#	+print() = Print a VTK file / Matplotlib visualization of the mesh (points and connections between nodes). Also print volumes.
class Mesh (object):
    def __init__(self):
        setDomain()

    def setDomain(self):
        pass

    def getPosition(self, ind):
        pass

    def getIndex(self, pos):
        pass

    def arrayToIndex(self, array):
        pass

    def indexToArray(self, ind):
        pass

    def vtkOrdering(self, array):
        pass

    def print(self):
        pass


#Mesh_2D_rm (Inherits from Mesh):
#
#Definition = Mesh class for a 2D rectangular mesh. The organization of the points will work as 0<=i<nx and 0<=j<ny. Also, for k parameter 0<=k<nPoints, k = nx*j+i.
#Attributes:
#	+xmin (double) = Left limit of the domain (closest to the Sun).
#	+xmax (double) = Right limit of the domain (farthest from the Sun).
#	+ymin (double) = Bottom limit of the domain.
#	+ymax (double) = Top limit of the domain.
#	+depth (double) = Artificial thickness of the domain, to make it three-dimensional.
#	+nx (int) = Number of nodes in the x direction.
#	+ny (int) = Number of nodes in the y direction.
#       +dx (float32) = Distance between adyacent horizontal nodes
#       +dy (float32) = Distance between adyacent vertical nodes
#       +boundaries ([Boundary]) = It is [Outer_2D_Rectangular].
#       +Mesh class attributes.
#Methods:
#	+Implementation of Mesh methods.
class Mesh_2D_rm (Mesh):
    def __init__(self, xmin, xmax, ymin, ymax, nx, ny, depth):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.nx = nx
        self.ny = ny
        self.depth = depth

        self.boundaries = [ob.Outer_2D_Rectangular(xmin, xmax, ymin, ymax)]
        self.xmin = self.boundaries[0].xmin
        self.xmax = self.boundaries[0].xmax
        self.ymin = self.boundaries[0].ymin
        self.ymax = self.boundaries[0].ymax
        self.dx = (self.xmax-self.xmin)/(self.nx-1)
        self.dy = (self.ymax-self.ymin)/(self.ny-1)
        self.setDomain()

#       +setDomain() = This function, with the values provided by the boundary files, will create the mesh, by setting up volumes, nPoints and any other subclass variable.
    def setDomain(self):
        self.nPoints = numpy.uint16(self.nx*self.ny)

        self.volumes = (self.dx*self.dy*self.depth)*numpy.ones((self.nPoints), dtype = 'float32')
        self.volumes[:self.nx] /= 2
        self.volumes[self.nx*(self.ny-1):] /= 2
        self.volumes[self.nx-1::self.nx] /= 2
        self.volumes[:self.nx*self.ny:self.nx] /= 2

        self.boundaries[0].location.extend(range(1, self.nx-1))
        self.boundaries[0].location.extend(range(0, self.nx*self.ny, self.nx))
        self.boundaries[0].location.extend(range(self.nx-1, self.nx*self.ny, self.nx))
        self.boundaries[0].location.extend(range(self.nx*(self.ny-1)+1, self.nx*self.ny-1))
        self.boundaries[0].location.sort()
        self.boundaries[0].location = numpy.asarray(self.boundaries[0].location, dtype = 'uint16')

#	+getPosition([int] i): [double, double y] = For a each index return its real position.
    def getPosition(self, ind):
        index2D = self.arrayToIndex(ind)
        return numpy.append(self.xmin+self.dx*index2D[:,0][:,None], self.ymin+self.dy*index2D[:,1][:,None], axis = 1)

#	+getIndex([double,double] pos): [double,double] = For each real position returns its index value.
    def getIndex(self, pos):
        indexes = numpy.zeros((numpy.shape(pos)[0],2))
        indexes[:,0] = (pos[:,0]-self.xmin)/self.dx
        indexes[:,1] = (pos[:,1]-self.ymin)/self.dy
        return indexes

#	+arrayToIndex([ind] array): [int, int] = For the indexes in the 1D array, obtain the indexes used for the particular mesh.
    def arrayToIndex(self, array):
        j, i = numpy.divmod(array, self.ny)
        return numpy.append(i[:,None], j[:,None], axis = 1)

#	+indexToArray([ind, ind] index): [int] = For the indexes used for the particular mesh, obtain the 1D version for the array.
    def indexToArray(self, ind):
        return ind[:,1]*self.nx+ind[:,0]

#       +vtkOrdering(array): array = The array received as argument is ordered in such a way it can be stored ina VTK file.
#           The result is returned as a new array.
    def vtkOrdering(self, array):
        dims = numpy.shape(array)
        if len(dims) == 1:
            return array.reshape((self.nx, self.ny, 1), order = 'F')
        else:
            tpl = tuple(numpy.reshape(copy.copy(array[:,i]),(self.nx, self.ny, 1), order = 'F') for i in range(dims[1]))
            if len(tpl) < 3:
                for i in range (3-len(tpl)):
                    tpl += (numpy.zeros_like(array[:,0].reshape((self.nx,self.ny,1))),)
            return tpl

#       +vtkOrdering(array): array = The array received as argument comes with vtk ordering and is reshaped to be stored properly in the code.
    def reverseVTKOrdering(self, array):
        dims = numpy.shape(array)
        if len(dims) == 1:
            return array.reshape((self.nPoints), order = 'F')
        else:
            return array.reshape((self.nPoints, 3), order = 'F')[:,:2]

#       +vtkReader(): Reader = Return the reader from module vtk that is able to read the vtk file.
    def vtkReader(self):
        return vtk.vtkXMLRectilinearGridReader()

#       +saveVTK(string filename, dict dictionary) = It calls the appropiate method in 'vtk' module to print the information of the system in a '.vtk' file.
    def saveVTK(self, filename, dictionary):
        i = numpy.arange(0, self.nx, dtype ='int16')
        j = numpy.arange(0, self.ny, dtype ='int16')
        ind = numpy.arange(0, self.nPoints, dtype = 'uint16')
        temp = numpy.zeros((1), dtype = 'int16')

        vtk.gridToVTK(filename, i, j, temp, pointData = dictionary)

#	+print() = Print a VTK file / Matplotlib visualization of the mesh (points and connections between nodes). Also print volumes.
    def print(self):
        i = numpy.arange(0, self.nx, dtype ='int16')
        j = numpy.arange(0, self.ny, dtype ='int16')
        ind = numpy.arange(0, self.nPoints, dtype = 'uint16')
        pos = self.getPosition(ind)
        temp = numpy.zeros((1), dtype = 'int16')

        cwd = os.path.split(os.getcwd())[0]
        vtkstring = cwd+'/results/mesh'
        vtk.gridToVTK(vtkstring, i, j, temp, pointData = {'volumes': self.vtkOrdering(self.volumes),\
                'positions': self.vtkOrdering(pos)})
                #'positions': (numpy.reshape(copy.copy(pos[:,0]),(self.nx,self.ny,1), order = 'F'), numpy.reshape(copy.copy(pos[:,1]),(self.nx,self.ny,1), order = 'F'), numpy.zeros((self.nx,self.ny,1)))})

#Data structures that define the information of the species in general
from vtk.util.numpy_support import vtk_to_numpy
import constants as c
import numpy

#Species (Abstract):
#
#Definition = This abstract class is the wrapper for everything related to particles. Any species of particles inherits from this class.
#Attributes:
#	+name (string) = some string descriptor that indicate the type/source of particles.
#       +dt (double) = timestep for the motion of the particle.
#	+q (double) = charge of the species.
#	+m (double) = mass of the species.
#	+q_over_m (double) = q/m.
#	+debye (double) = Debye length of the spcies.
#	+spwt (int) = specific weight of the species.
#       +pos_dim (int) = number of components of position
#       +vel_dim(int) = number of components of the velocity
#	+mesh_values (Particles_In_Mesh) = object to store anything related with the particles in the mesh.
#	+part_values (Particles) = oject to store anything related to the actual particles in physical space.
#Notes:
#       +__init__() receives nPoints from mesh
#       +saveVTK(Mesh mesh): dictionary = Return the attributes of the Species to be printed in the VTK file.
#           The process is handled through Particles_In_Mesh.
#       +loadVTK(Mesh mesh, output) = Takes information of the species from a VTK file through 'output' and stores it in the corresponding attributes.
#           The process is handled through Particles_In_Mesh.
class Species(object):
    def __init__(self, name, dt, n_q, n_m, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints, n_num_tracked = 0):
        self.name = name
        self.dt = dt
        self.q = n_q
        self.m = n_m
        self.q_over_m = self.q/self.m
        self.debye = n_debye
        self.spwt = n_spwt
        self.pos_dim = n_pos_dim
        self.vel_dim = n_vel_dim
        self.mesh_values = Particles_In_Mesh(n_nPoints, n_vel_dim)
        self.part_values = Particles(n_max_n, n_pos_dim, n_vel_dim, n_num_tracked)

    def saveVTK(self, mesh):
        return self.mesh_values.saveVTK(mesh, self.name)

    def loadVTK(self, mesh, output):
        self.mesh_values.loadVTK(mesh, output, self.name)
        mesh.loadSpeciesVTK(self)
    
    def saveParticlesTXT(self):
        data, names, array = self.part_values.saveParticlesTXT()
        species_data = self.name+"\t"+data
        return  species_data, names, array


#Particles_In_Mesh (Abstract)(Composition with Species):
#
#Definition = Store values of distributions in the mesh related with particles. It wraps all the attributes that all particles must have.
#Attributes:
#	+nPoints (int) = Number of nodes in the mesh (Same as mesh class).
#	+density ([double]) = Density values at each node.
#	+velocity ([double, double]) = Velocity at each node. Rows are different points, columns are (x,y,z) components if they are available.
#	+temperature ([double]) = Temperature values at each node.
#       +residuals([double]) = remnants from injection of particles at the previous step.
#Methods:
#       +saveVTK(Mesh mesh, string name): dictionary = Return the attributes of the Species to be printed in the VTK file.
#       +loadVTK(Mesh mesh, output, string name) = Takes information of the species from a VTK file through 'output' and stores it in the corresponding attributes.
class Particles_In_Mesh(object):
    def __init__(self, n_nPoints, n_vel_dim):
        self.nPoints = n_nPoints
        self.density = numpy.zeros((self.nPoints))
        self.velocity = numpy.zeros((self.nPoints, n_vel_dim))
        self.temperature = numpy.zeros((self.nPoints))
        self.residuals = numpy.zeros((self.nPoints))

    def saveVTK(self, mesh, name):
        return {name+"-density" : mesh.vtkOrdering(self.density),\
                name+"-velocity": mesh.vtkOrdering(self.velocity),\
                name+"-temperature": mesh.vtkOrdering(self.temperature/c.EV_TO_K),\
                name+"-residuals": mesh.vtkOrdering(self.residuals)}

    def loadVTK(self, mesh, output, name):
        self.density = mesh.reverseVTKOrdering(vtk_to_numpy(output.GetPointData().GetArray(name+"-density")))
        self.velocity = mesh.reverseVTKOrdering(vtk_to_numpy(output.GetPointData().GetArray(name+"-velocity")))
        self.temperature = mesh.reverseVTKOrdering(vtk_to_numpy(output.GetPointData().GetArray(name+"-temperature")))*c.EV_TO_K


#Particles(Composition with Species):
#
#Definition = Stores values related with the particles themselves.
#Attributes:
#	+current_n (int) = current number of particles.
#	+max_n (int) = max. number of particles for the species.
#	+position ([double,double]) = Position of every particle. Rows are different particles, columns are (x,y,z) if available.
#	+velocity ([double,double]) = Position of every particle. Rows are different particles, columns are (x,y,z) components if available.
#	+num_tracked (int) = Size of particles being tracked. Defaults to 0 meaning that the species is not being tracked.
#	+trackers ([int]) = array of size num_tracked that store the indices of the particles as stored in positions.
class Particles(object):
    def __init__(self, n_max_n, n_pos_dim, n_vel_dim, num_tracked):
        self.current_n = numpy.uint32(0)
        self.max_n = numpy.uint32(n_max_n)
        self.position = numpy.zeros((n_max_n, n_pos_dim))
        self.velocity = numpy.zeros((n_max_n, n_vel_dim))
        self.num_tracked = numpy.uint16(num_tracked)
        if num_tracked != 0:
            self.trackers = self.max_n*numpy.ones((num_tracked), dtype = numpy.uint32)

    def saveParticlesTXT(self):
        names = "position\tvelocity"
        data = "{:6d}-{:1d}-{:1d}".format(self.current_n, numpy.shape(self.position)[1], numpy.shape(self.velocity)[1])
        array = numpy.append(self.position[:self.current_n,:], self.velocity[:self.current_n,:], axis = 1)
        return data, names, array

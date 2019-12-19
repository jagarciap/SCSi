#Data structures that define the information of the species in general
import numpy

#Species (Abstract):
#
#Definition = This abstract class is the wrapper for everything related to particles. Any species of particles inherits from this class.
#Attributes:
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
#       +__init__() also receives nPoints from mesh
class Species(object):
    def __init__(self, n_q, n_m, n_q_over_m, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints):
        self.q = nq_q
        self.m = n_m
        self.q_over_m = n_q_over_m
        self.debye = n_debye
        self.spwt = n_spwt
        self.pos_dim = n_pos_dim
        self.vel_dim = n_vel_dim
        self.mesh_values = Particles_In_Mesh(n_nPoints, n_vel_dim)
        self.part_values = Particles(n_max_n, n_pos_dim, n_vel_dim)


#Particles_In_Mesh (Abstract)(Composition with Species):
#
#Definition = Store values of distributions in the mesh related with particles. It wraps all the attributes that all particles must have.
#Attributes:
#	+nPoints (int) = Number of nodes in the mesh (Same as mesh class).
#	+density ([double]) = Density values at each node.
#	+velocity ([double, double]) = Velocity at each node. Rows are different points, columns are (x,y,z) components if they are available.
#       +residuals([double]) = remnants from injection of particles at the previous step.
class Particles_In_Mesh(object):
    def __init__(self, n_nPoints, n_vel_dim):
        self.nPoints = n_nPoints
        self.density = numpy.zeros((nPoints))
        self.velocity = numpy.zeros((nPoints, n_vel_dim))
        self.residuals = numpy.zeros((nPoints))


#Particles:
#
#Definition = Stores values related with the particles themselves.
#Attributes:
#	+current_n (int) = current number of particles.
#	+max_n (int) = max. number of particles for the species.
#	+position ([double,double]) = Position of every particle. Rows are different particles, columns are (x,y,z) if available.
#	+velocity ([double,double]) = Position of every particle. Rows are different particles, columns are (x,y,z) components if available.
class Particles(object):
    def __init__(self, n_max_n, n_pos_dim, n_vel_dim):
        self.current_n = numpy.uint32(0)
        self.max_n = n_max_n
        self.position = numpy.zeros((n_max_n, n_pos_dim))
        self.velocity = numpy.zeros((n_max_n, n_vel_dim))

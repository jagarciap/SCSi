#Data structures to hold domain information
import numpy
import constants as c

#Mesh (Abstract)(Association between Mesh and PIC):
#
#Definition = Defines the type of mesh.
#Attributes:
#	+nPoints (int) = Number of points in the mesh.
#	+volumes ([double]) = Volume of each node.
#	+pic (PIC) = Object that provides pic methods (Association between Mesh and PIC).
#Methods:
#       +setDomain() = This function, with the values provided by the boundary files, will create the mesh, by setting up volumes, nPoints and any other subclass variable.
#	+getPosition([int] i): [double, double y] = For a each index return its real position.
#	+getVolumes() = calculates the volume of each node and stores it in volumes attribute.
#	+print() = Print a VTK file / Matplotlib visualization of the mesh (points and connections between particles).
class Mesh (object):
    def __init__(self, pic_object):
        self.pic = pic_object
        setDomain()

    def setDomain(self):
        pass

    def getPosition(self, ind):
        pass

    def getVolumes(self):
        pass

    def print(self):
        pass


#Mesh_2D_rm (Inherits from Mesh):
#
#Definition = Mesh class for a 2D rectangular mesh.
#Attributes:
#	+xmin (double) = Left limit of the domain (closest to the Sun).
#	+xmax (double) = Right limit of the domain (farthest from the Sun).
#	+ymin (double) = Bottom limit of the domain.
#	+ymax (double) = Top limit of the domain.
#	+depth (double) = Artificial thickness of the domain, to make it three-dimensional.
#	+nx (int) = Number of nodes in the x direction.
#	+ny (int) = Number of nodes in the y direction.
#Methods:
#	+Implementation of Mesh methods.
class Mesh_2D_rm (Mesh):
    def __init__(self, pic_object):
        self.


class Domain(object):
    def __init__(self,nz,nr,dz,dr,zmin,rmin,dt):
        self.nz=nz				#number of nodes*/
        self.nr=nr
        self.dz=dz
        self.dr=dr			#cell spacing*/
        self.zmin = zmin
        self.rmin = rmin
        self.dt=dt
        
        self.phi=numpy.zeros((nz,nr))		#potential*/
        self.ef =numpy.zeros((2,nz,nr))		#electric field*/

        self.i_vel = numpy.zeros((2,nz,nr))	#ion velocity
        self.i_current = numpy.zeros((nz,nr))   #current created by ions
        self.e_vel = numpy.zeros((2,nz,nr))
        self.e_current = numpy.zeros((nz,nr))   #current created by electrons

        self.ndi=numpy.zeros((nz,nr))		#ion density*/
        self.rho=numpy.zeros((nz,nr))		#charge density*/
        self.nde=numpy.zeros((nz,nr))           #neutrals density

        self.node_type=numpy.zeros((nz,nr), dtype = numpy.int8)		#node type for geometry
        self.node_volume= dz*dr*numpy.ones((nz,nr))   #node volumes

        self.temperature = numpy.zeros((nz,nr))         #Temperature in each node
        self.mf = numpy.zeros((nz,nr))
        self.mu_e = numpy.zeros((nz,nr))
        
class Particle(object):
    def __init__(self,num):
        self.x=numpy.zeros((num,2)) #Position
        self.v=numpy.zeros((num,2)) #Velocity
        self.printable = numpy.zeros(num, dtype = numpy.int8) #Print in Particle tracker
        self.print_ind = [] #storage of printing indexes
		
class Species(object):
    def __init__(self):
        self.mass=[]
        self.charge=[]
        self.charge_over_mass=[]
        self.spwt=[]

        self.np=[]
        self.np_alloc=[]
        self.part=[]

class NeumannBoundaries(object):
    def __init__(self,world):
        world.node_type[c.Z_IND:,0] = 5
        world.node_type[-1,:] = 5
        world.node_type[c.Z_IND:,-1] = 5

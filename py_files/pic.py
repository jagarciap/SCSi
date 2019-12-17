#Data structures that contain PIC numerical methods
import numpy
import constants as c

#PIC (Abstract)(Association between Mesh and PIC):
#
#Definition = Indicate the methods that all PIC classes have to implement. Each PIC concrete class will depend on the type of mesh, as well as the type of PIC algorithm implemented.
#Methods:
#	+scatter([double, double] positions, [double] values, [double] field) = Receives the positions of the particles, and makes scatter procedure, calculating the values of field for each node.
#	+gather([double, double] positions, [double] field): [double]field_p = Calculates values of the field in particles' positions, returning these values in an array as long as positions, called field_p.
class PIC (object):
    def scatter(self, positions, values, field):
        pass

    def gather(self, positions, field):
        pass


#PIC_2D_rm1o (Inherits from PIC):
#
#Definition = PIC class for rm10 ('r'ectangular 'm'esh, '1'st 'o'rder implementation).
#Methods:
#	+Implementation of PIC class methods.
#	+scatter_density = return densities of that species in every node of the mesh.
#	+scatter_velocity = return velocities of that species in every node of the mesh.
#	+scatter_flux = return flux of particles of that species into every indicated node (not all the mesh).
class PIC_2D_rm1o(PIC):
    string = "PIC method for a rectangular 2D mesh, making interpolations at first order"


def scatter_opt(lc, value, field):
    lc = lc[numpy.logical_and(lc[:,0]>0,lc[:,1]>0),:]
    i = lc[:,0].astype(int)
    j = lc[:,1].astype(int)
    di = lc[:,0] - i
    dj = lc[:,1] - j
    
    nx, ny = numpy.shape(field)
    base = numpy.logical_and(i < nx, j < ny)
    ind0 = numpy.flatnonzero(base)
    indx = numpy.flatnonzero(numpy.logical_and(i+1 < nx, base))
    indy = numpy.flatnonzero(numpy.logical_and(j+1 < ny,base))
    indxy = numpy.flatnonzero(numpy.logical_and(i+1 < nx, j+1 < ny))

    numpy.add.at(field, (i[ind0],j[ind0]),(1-di[ind0])*(1-dj[ind0])*value[ind0])
    numpy.add.at(field, (i[indx]+1,j[indx]), di[indx]*(1-dj[indx])*value[indx])
    numpy.add.at(field, (i[indy], j[indy]+1), (1-di[indy])*dj[indy]*value[indy])
    numpy.add.at(field, (i[indxy]+1, j[indxy]+1), di[indxy]*dj[indxy]*value[indxy])

def ScatterSpecies (world,species,field):
    zmin = world.zmin
    rmin = world.rmin
    dz = world.dz
    dr = world.dr
    np = species.np

    data = species.part.x[:np,:]
    
    lc = numpy.zeros_like(data)
    lc[:,0] = (data[:,0]-zmin)/dz
    lc[:,1] = (data[:,1]-rmin)/dr


    field*=0
   
    #scatter particles to the mesh
    scatter_opt(lc, species.spwt*numpy.ones((np)),field)
    
    #divide by cell volume
    field/=world.node_volume
    print(numpy.max(field))

def ScatterSpeed (world,species,field):
    zmin = world.zmin
    rmin = world.rmin
    dz = world.dz
    dr = world.dr
    np = species.np

    data = species.part.x[:np,:]
    
    lc = numpy.zeros_like(data)
    lc[:,0] = (data[:,0]-zmin)/dz
    lc[:,1] = (data[:,1]-rmin)/dr

    field*=0
    
    #scatter particles to the mesh
    for dim in range(2):
        scatter_opt(lc, species.part.v[:np,dim],field[dim])
    
    if species.charge == -c.QE:
        field*=numpy.where(world.ndi==0,0,species.spwt/world.ndi/world.node_volume)
    elif species.charge == c.QE:
        field*=numpy.where(world.nde==0,0,species.spwt/world.nde/world.node_volume)
    
def gather_opt(data,lc):
    i = lc[:,0].astype(int)
    j = lc[:,1].astype(int)
    di = lc[:,0] - i
    dj = lc[:,1] - j

    n = numpy.shape(i)[0]
    nx, ny = numpy.shape(data)
    base = numpy.logical_and(i < nx, j < ny)
    ind0 = numpy.flatnonzero(base)
    indx = numpy.flatnonzero(numpy.logical_and(i+1 < nx, base))
    indy = numpy.flatnonzero(numpy.logical_and(j+1 < ny,base))
    indxy = numpy.flatnonzero(numpy.logical_and(i+1 < nx, j+1 < ny))

    a = numpy.zeros((n))
    numpy.add.at(a, ind0, data[i[ind0],j[ind0]]*(1-di[ind0])*(1-dj[ind0]))
    numpy.add.at(a, indx, data[i[indx]+1,j[indx]]*di[indx]*(1-dj[indx]))
    numpy.add.at(a, indy, data[i[indy],j[indy]+1]*(1-di[indy])*dj[indy])
    numpy.add.at(a, indxy, data[i[indxy]+1,j[indxy]+1]*di[indxy]*dj[indxy])

    return a

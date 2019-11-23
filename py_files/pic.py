import numpy
import constants as c

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

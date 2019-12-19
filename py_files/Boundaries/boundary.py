#Data structure of the Boundaries
import numpy

#Boundary (Abstract-like: will also have common methods that can be used by sub-classes. In composition with mesh):
#
#Definition = Class that shows the methods and attributes needed for each particular boundary arrangement.
#Attributes:
#	+type (string) = string indicating the type of boundary. It will be of the form "[Inner or Outer] - [Sorce, e.g. Spacecraft, Component, etc.]".
#	+location ([int]) = array containing the indices of the nodes that the boundary represents.
#Methods:
#	+applyElectricBoundary(Electric_Field) = Applies the boundary condition to the electric field passed as argument.
#	+applyMagneticBoundary(Magnetic_Field) = Applies the boundary condition to the magnetic field passed as argument.
#	+applyParticleBoundary(Species) = Applies the boundary condition to the species passed as argument.
class Boundary(object):
    def applyElectricBoundary(self, e_field):
        pass
    
    def applyMagneticBoundary(self, m_field):
        pass

    def applyParticleBoundary(self, species):
        pass

#       +Add new particles to species. pos and vel represent their positions and velocities, respectively.
    def addParticles (self, species, pos, vel):
        n = numpy.shape(pos)[0]
        if species.part_values.current_n + n > species.part_values.max_n:
            raise ValueError ("Too many particles")
    
        #store position and velocity of this particle
        species.part.x[species.part_values.current_n:species.part_values.current_n+n]= pos
        species.part.v[species.part_values.current_n:species.part_values.current_n+n]= vel
        
        #increment particle counter
        species.part_values.current_n += n

#       +Function that inject particles into the domain.
#       +Parameters: 
#       +location ([ind]) = Nodes indicating the faces where the particles are going to be inserted. Each node represents the volume surrounded it. Location should be ordered increasingly.
#       +pic (PIC) = Instance of PIC for calculations. It also constains mesh, which is used in the function.
#       +species (Species) = Species to be inserted. It contains inside Particles_In_Mesh 'residuals' which is for each node the values of remnants from the insertion at the previous step.
#       +delta_n ([double]) = For each node in location, the density that is going to be inserted at this timestep. The array is ordered with the order of the nodes in mesh.
#       +n_vel ([double,double]) = For each node in location, the thermal velocity (for a the MB distribution) that the inserted particles will represent. Same order as delta_n.
def injectParticlesFace(self, location, pic, species, delta_n, n_vel): 
    # Floating point production rate
    mpf_new = delta_n*pic.mesh.volumes[location]/species.spwt+species.mesh_values.residuals[location]
    # Truncate down, adding randomness
    mp_new = (mpf_new+numpy.random.rand(len(location))).astype(int)
    # Save fraction part
    species.mesh_values.residuals[location] = mpf_new-mp_new
    # Generate this many particles
    total_new = numpy.sum(mp_new)
    pos = numpy.zeros((total_new, species.pos_dim))
    vel = numpy.zeros((total_new, species.vel_dim))
    #Indices to operate over
    ind = []
    ind.extend(for i in lo i*numpy.ones((mp_new)))
    ind = 



    i = c.NR1
    count = 0
    while i < c.NR2:
    
        lc=numpy.zeros(2)
        lc[1]=i+0.5
        
        dni=dn*Gather(world.node_volume,lc)             
        # This exists to fix the problem with node volumes on the edges
        if i == 0:
            dni += dn*0.5*world.node_volume[0,i]
        elif i == world.nr-2:
            dni += dn*0.5*world.node_volume[0,i]
        
        #floating point production rate
        mpf_new = dni/species.spwt + mpf_rem[i]
        
        #truncate down, adding randomness
        mp_new = int(numpy.trunc(mpf_new+rand.random()))
        
        #save fraction part
        mpf_rem[i] = mpf_new - mp_new    #new fractional reminder
        
        #generate this many particles
        pos=numpy.zeros((mp_new,2))
        # MODIFIED
        pos[:,0] += numpy.random.rand(mp_new)*(nVel*world.dt)
        pos[:,1] = world.rmin+(i+numpy.random.rand(mp_new))*world.dr
        
        #Determine velocity of neutrals
        vel = sampleIsotropicVel(nVel, mp_new)
        vel[:,0] = numpy.abs(vel[:,0])
                   
        AddParticle(species,pos,vel)
        count +=mp_new
        i +=1

    print("Injected particles: ",count)
    return count

def sampleIsotropicVel(vth, num):
    #pick a random angle
    theta = 2*numpy.pi*numpy.random.rand(num)
    
    #pick a random direction for n[2]
    n = numpy.zeros((num,2))
    n[:,0] = numpy.cos(theta)
    n[:,1] = numpy.sin(theta)
    
    #pick maxwellian velocities
    rand_spread = numpy.random.rand(num,3)
    vm = numpy.sqrt(2)*vth*(rand_spread[:,0]+rand_spread[:,1]+rand_spread[:,2]-1.5)
    
    vel = n*vm[:,None]
    return vel

    # Particles going out of domain 
    ind = numpy.flatnonzero(numpy.logical_or(numpy.logical_or(numpy.logical_or(species.part.x[:np,0]< zmin, \
            species.part.x[:np,0] > (((nz-1)*dz)+zmin)), \
            species.part.x[:np,1] > rmax), \
            species.part.x[:np,1] < rmin))
    count2 = numpy.shape(ind)[0]
    ndi_out = RemoveParticle(species,ind)

def RemoveParticle(species,ind):
    dz = c.DZ
    num = numpy.shape(numpy.where(species.part.x[ind,0] > dz))[1]

    temp = species.np
    species.np -= numpy.shape(ind)[0]
    species.part.x[:species.np,:] = numpy.delete(species.part.x[:temp,:], ind, axis = 0)
    species.part.v[:species.np,:] = numpy.delete(species.part.v[:temp,:], ind, axis = 0)
    return num


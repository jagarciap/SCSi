#Data structure of the Boundaries
import numpy
import pdb

from Species.species import Species

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
#       +createDummyBox(self, location, pic, species, delta_n, n_vel, shift_vel) = create the dummy box with particles in it.
#       +addParticles(Species species, [double, double] pos, [double, double] vel) = Add to Species the new particles, each represented by a row in pos and vel.
#       +removeParticles(Species species, [ind] ind, Boolean tracker) = Removes the particles from species stored at 'ind' positions. tracker indicates whether there is need for handling a Tracker instance.
#       +sampleIsotropicVelocity(double vth, int num) = It receives the most probable speed vth = \sqrt{2kT/m} and creates num random 2D velocities with their magnitudes following a Maxwellian distribution.
class Boundary(object):
    def applyElectricBoundary(self, e_field):
        pass
    
    def applyMagneticBoundary(self, m_field):
        pass

    def applyParticleBoundary(self, species):
        pass

    def createDummyBox (self, location, pic, species, delta_n, n_vel, shift_vel):
        pass

#       +It receives the most probable speed vth = \sqrt{2kT/m} and creates num random 2D velocities with their magnitudes following a Maxwellian distribution.
#       NOTE: This function needs to be revised. random should not spread throughout different cells, and should use the same seed for every function call.
#    @nb.vectorize(signature = nb.double[:], target='cpu')
    def sampleIsotropicVelocity(self, vth, num):
        #pick a random angle
        theta = 2*numpy.pi*numpy.random.rand(num)
        #pick a random direction for n[2]
        n = numpy.append(numpy.cos(theta)[:,None], numpy.sin(theta)[:,None], axis = 1)
        #pick maxwellian velocities
        rand_spread = numpy.random.rand(num,3)
        vm = numpy.sqrt(2)*vth*(rand_spread[:,0]+rand_spread[:,1]+rand_spread[:,2]-1.5)
        #2D components of velocity 
        vel = n*vm[:,None]
        return vel

#       +Add new particles to species. pos and vel represent their positions and velocities, respectively.
    def addParticles(self, species, pos, vel):
        n = numpy.shape(pos)[0]
        if species.part_values.current_n + n > species.part_values.max_n:
            raise ValueError ("Too many particles")
        #store position and velocity of this particle
        species.part_values.position[species.part_values.current_n:species.part_values.current_n+n]= pos
        species.part_values.velocity[species.part_values.current_n:species.part_values.current_n+n]= vel
        #increment particle counter
        species.part_values.current_n += n

#       +Eliminates the particles  denoted with the indices ind. 
    def removeParticles(self, species, ind):
        #Eliminating particles
        temp = species.part_values.current_n
        species.part_values.current_n -= numpy.shape(ind)[0]
        species.part_values.position[:species.part_values.current_n,:] = numpy.delete(species.part_values.position[:temp,:], ind, axis = 0)
        species.part_values.velocity[:species.part_values.current_n,:] = numpy.delete(species.part_values.velocity[:temp,:], ind, axis = 0)
        #Updating tracker
        #NOTE: Include unsorted to sorted arrays
        if species.part_values.num_tracked != 0:
            ind_c = 0
            tracker_c = 0
            n = 0
            while ind_c != len(ind) and tracker_c != species.part_values.num_tracked:
                if ind1[ind_c] < trackers1[tracker_c]:
                    n += 1
                    ind_c += 1
                    break
                elif ind1[ind_c] > trackers1[tracker_c]:
                    trackers1[tracker_c] -= n
                    tracker_c += 1
                    break
                elif ind1[ind_c] == trackers1[tracker_c]:
                    trackers1[tracker_c] = species.part_values.max_n
                    n += 1
                    tracker_c += 1
                    ind_c += 1
            if ind_c == len(ind1) and tracker_c < len(trackers1):
                trackers1[tracker_c:] -= n

#       +Function that inject particles into the domain.
#       +Parameters: 
#       +location ([ind]) = Nodes indicating the faces where the particles are going to be inserted. Each node represents the volume surrounding it. Location should be ordered increasingly.
#       +pic (PIC) = Instance of PIC for calculations. It also constains mesh, which is used in the function.
#       +species (Species) = Species to be inserted. It contains inside Particles_In_Mesh 'residuals' which is for each node the values of remnants from the insertion at the previous step.
#       +delta_n ([double]) = For each node in location, the density that is going to be inserted at this timestep. The array is ordered with the order of the nodes in mesh.
#       +n_vel ([double,double]) = For each node in location, the thermal velocity (for a the MB distribution) that the inserted particles will represent. Same order as delta_n.
#       +shift_vel ([double, double]) = For each node in location, an added velocity that does not come from temperature origin (like solar wind average speed.
#NOTE: I am doing a huge change in the conception of inject particle. I should give this a second thinking
    def injectParticlesFace(self, location, pic, species, delta_n, n_vel, shift_vel): 
        # Useful constant
        box_x = shift_vel[:,0]*species.dt
        # Floating point production rate
        mpf_new = delta_n*pic.mesh.volumes[location]*box_x/pic.mesh.dx/species.spwt+species.mesh_values.residuals[location]
        # Truncate down, adding randomness
        mp_new = (mpf_new+numpy.random.rand(len(location))).astype(int)
        # Save fraction part
        species.mesh_values.residuals[location] = mpf_new-mp_new
        # Generate this many particles
        total_new = numpy.sum(mp_new)
        pos = numpy.zeros((total_new, species.pos_dim))
        vel = numpy.zeros((total_new, species.vel_dim))
        # Preparing positions (this part is mesh-dependant [To be depurated later or changed each time].
        # Preparating velocities. This part is mesh-dependant trough sampleIsotropicVelocity.
        phys_loc = pic.mesh.getPosition(location)
        c = 0
        for i in range(len(location)):
            vel[c:c+mp_new[i],:] += self.sampleIsotropicVelocity(n_vel[i], mp_new[i])+shift_vel[i,:]
            pos[c:c+mp_new[i],0] += phys_loc[i,0] + numpy.random.rand(mp_new[i])*box_x[i]
            #pos[c:c+mp_new[i],0] += phys_loc[i,0] + (numpy.random.rand(mp_new[i])-0.5)*pic.mesh.dx
            pos[c:c+mp_new[i],1] += phys_loc[i,1] + (numpy.random.rand(mp_new[i])-0.5)*pic.mesh.dy
            c += mp_new[i]

        # Adjusting boundaries
        pos[:,1] += numpy.where(pos[:,1] <= pic.mesh.ymin, pic.mesh.dy/1.99, 0.0)
        pos[:,1] -= numpy.where(pos[:,1] >= pic.mesh.ymax, pic.mesh.dy/1.99, 0.0)
        pos[:,0] += numpy.where(pos[:,0] <= pic.mesh.xmin, pic.mesh.dx*0.001, 0.0)
        #pos[:,0] -= numpy.where(pos[:,0] >= pic.mesh.xmax, pic.mesh.dx/1.99, 0.0)
        #vel[:,0] = numpy.abs(vel[:,0])
        #Adding particles
        self.addParticles(species,pos,vel)
        print("Injected particles: ",total_new)
        print("Total{}".format(species.type), species.part_values.current_n)

#       +injectParticlesDummyBox([int] location, PIC pic, Species species, [double] delta_n, [double] n_vel, [double] shift_vel) = Inject the particles in location indices by creating dummy boxes around them, creating particles
#       	inside of them, moving the particles, and then adding the ones that entered into the computational domain.
    def injectParticlesDummyBox(self, location, pic, species, delta_n, n_vel, shift_vel):
        # Creating temporary species
        ghost = Species(species.dt, species.q, species.m, species.debye, species.spwt, int(species.part_values.max_n/10), species.pos_dim, species.vel_dim, species.mesh_values.nPoints)
        self.createDummyBox(location, pic, ghost, delta_n, n_vel, shift_vel)
        np = ghost.part_values.max_n
        ghost.part_values.position[:np,:] += ghost.part_values.velocity[:np,:]*ghost.dt
        ind = numpy.flatnonzero(numpy.logical_and(ghost.part_values.position[:np,0] > pic.mesh.xmin,\
                                numpy.logical_and(ghost.part_values.position[:np,0] < pic.mesh.xmax,\
                                numpy.logical_and(ghost.part_values.position[:np,1] > pic.mesh.ymin, ghost.part_values.position[:np,1] < pic.mesh.ymax))))
        self.addParticles(species, ghost.part_values.position[ind,:], ghost.part_values.velocity[ind])
        print("Injected particles: ",len(ind))
        print("Total{}".format(species.type), species.part_values.current_n)

import constants as c
import numpy
import pdb

from Boundaries.boundary import Boundary

#Outer_2D_Rectangular (Inherits from Boundary):
#
#Definition = Outer boundary for a rectangular mesh
#Attributes:
#	+type (string) = "Outer - Domain"
#	+xmin (double) = Left limit of the domain (closest to the Sun).
#	+xmax (double) = Right limit of the domain (farthest from the Sun).
#	+ymin (double) = Bottom limit of the domain.
#	+ymax (double) = Top limit of the domain.
#	+location ([int]) = Boundary attribute.
#Methods:
#	+Boundary methods.
class Outer_2D_Rectangular(Boundary):
    def __init__(self, x_min, x_max , y_min, y_max):
        self.type = "Outer - Domain"
        self.xmin = x_min
        self.xmax = x_max
        self.ymin = y_min
        self.ymax = y_max
        self.location = []

##	+applyElectricBoundary(Electric_Field) = Applies the boundary condition to the electric field passed as argument. So far a 0V Dirichlet boundary condition is applied.
##NOTE: Modifified to add a constant -20V at the right boundary (2020_03_12)
#    def applyElectricBoundary(self, e_field):
#        #Location of dirichlet and neumann boundaries
#        dirichlet_loc_1 = numpy.arange(c.NX-1, c.NX*c.NY, c.NX)
#        neumann_loc = numpy.delete(self.location, numpy.arange(c.NX-1,c.NX+(c.NY-1)*2, 2))
#        dirichlet_loc_2 = numpy.arange(0, c.NX*c.NY, c.NX)
#        #Values
#        dirichlet_val_1 = -20*numpy.ones_like(dirichlet_loc_1)
#        dirichlet_val_2 = numpy.zeros_like(dirichlet_loc_2)
#        neumann_val = numpy.zeros_like(neumann_loc)
#        #Applying values
#        e_field.dirichlet(dirichlet_loc_1, dirichlet_val_1)
#        e_field.dirichlet(dirichlet_loc_2, dirichlet_val_2)
#        e_field.neumann(neumann_loc, neumann_val)

#	+applyElectricBoundary(Electric_Field) = Applies the boundary condition to the electric field passed as argument. So far a 0V Dirichlet boundary condition is applied.
    def applyElectricBoundary(self, e_field):
        values = numpy.zeros((len(self.location)))
        e_field.dirichlet(self.location, values)

#	+applyMagneticBoundary(Magnetic_Field) = Applies the boundary condition to the magnetic field passed as argument.
#       No magnetic field so far
    def applyMagneticBoundary(self, m_field):
        pass

#	+applyParticleBoundary(Species) = Applies the boundary condition to the species passed as argument.
#           type_boundary indicates the type of boundary method to apply to particles. 'open', the default mthod, deletes them. 'reflective' reflects them back to the dominion.
#           **kwargs may contain arguments necessary for inner methods.
    def applyParticleBoundary(self, species, type_boundary, **kwargs):
        if type_boundary == 'open':
            self.applyParticleOpenBoundary(species)
        elif type_boundary == 'reflective':
            self.applyParticleReflectiveBoundary(species, kwargs['old_species'])
        else:
            raise ValueError("Called invalid boundary method")

#       +applyParticleOpenBoundary(Species) = Deletes particles at or outside of the boundaries.
    def applyParticleOpenBoundary(self, species):
        #Just for convenience in writing
        np = species.part_values.current_n
        #Finding the particles
        ind = numpy.flatnonzero(numpy.logical_or(numpy.logical_or(numpy.logical_or(species.part_values.position[:np,0] <= self.xmin, \
                species.part_values.position[:np,0] >= self.xmax), \
                species.part_values.position[:np,1] >= self.ymax), \
                species.part_values.position[:np,1] <= self.ymin))
        # Eliminating particles
        super().removeParticles(species,ind)
        count2 = numpy.shape(ind)[0]
        print('Number of {} eliminated:'.format(species.name), count2)

#       +applyParticleReflectiveBoundary(Species species, Species old_species) = Reflects the particles back into the domain.
#           old_species refers to the state of species in the previous step.
    def applyParticleReflectiveBoundary(self, species, old_species):
        #For ease of typing
        np = species.part_values.current_n
        delta = 1e-5
        #Upper wall
        ind = numpy.flatnonzero(species.part_values.position[:np,1] >= self.ymax)
        species.part_values.velocity[ind,1] *= -1
        pos = (self.ymax-delta)*numpy.ones((len(ind), species.pos_dim))
        pos[:,0] = (species.part_values.position[ind,0]-old_species.part_values.position[ind,0])/(species.part_values.position[ind,1]-old_species.part_values.position[ind,1])\
                *(self.ymax-delta-old_species.part_values.position[ind,1])+old_species.part_values.position[ind,0]
        species.part_values.position[ind,:] = pos
        #Left wall
        ind = numpy.flatnonzero(species.part_values.position[:np,0] <= self.xmin)
        ind_v = numpy.flatnonzero(species.part_values.position[ind,1] == self.ymax-delta)
        species.part_values.velocity[ind[ind_v],1] *= numpy.where(species.part_values.position[ind[ind_v],0] == self.xmin+delta, 1, -1)
        species.part_values.velocity[ind,0] *= -1
        pos = (self.xmin+delta)*numpy.ones((len(ind), species.pos_dim))
        pos[:,1] = (species.part_values.position[ind,1]-old_species.part_values.position[ind,1])/(species.part_values.position[ind,0]-old_species.part_values.position[ind,0])\
                *(self.xmin+delta-old_species.part_values.position[ind,0])+old_species.part_values.position[ind,1]
        species.part_values.position[ind,:] = pos
        #Right wall
        ind = numpy.flatnonzero(species.part_values.position[:np,0] >= self.xmax)
        ind_v = numpy.flatnonzero(species.part_values.position[ind,1] == self.ymax-delta)
        species.part_values.velocity[ind[ind_v],1] *= numpy.where(species.part_values.position[ind[ind_v],0] == self.xmax-delta, 1, -1)
        species.part_values.velocity[ind,0] *= -1
        pos = (self.xmax-delta)*numpy.ones((len(ind), species.pos_dim))
        pos[:,1] = (species.part_values.position[ind,1]-old_species.part_values.position[ind,1])/(species.part_values.position[ind,0]-old_species.part_values.position[ind,0])\
                *(self.xmax-delta-old_species.part_values.position[ind,0])+old_species.part_values.position[ind,1]
        species.part_values.position[ind,:] = pos
        #Lower wall
        ind = numpy.flatnonzero(species.part_values.position[:np,1] <= self.ymin)
        ind_v = numpy.flatnonzero(species.part_values.position[ind,0] == self.xmax-delta)
        species.part_values.velocity[ind[ind_v],0] *= numpy.where(species.part_values.position[ind[ind_v],1] == self.ymin+delta, 1, -1)
        ind_v = numpy.flatnonzero(species.part_values.position[ind,0] == self.xmin+delta)
        species.part_values.velocity[ind[ind_v],0] *= numpy.where(species.part_values.position[ind[ind_v],1] == self.ymin+delta, 1, -1)
        species.part_values.velocity[ind,1] *= -1
        pos = (self.ymin+delta)*numpy.ones((len(ind), species.pos_dim))
        pos[:,0] = (species.part_values.position[ind,0]-old_species.part_values.position[ind,0])/(species.part_values.position[ind,1]-old_species.part_values.position[ind,1])\
                *(self.ymin+delta-old_species.part_values.position[ind,1])+old_species.part_values.position[ind,0]
        species.part_values.position[ind,:] = pos

#       +createDummyBox([ind]location, PIC pic, Species species, [double] delta_n, [double] n_vel, [double] shift_vel) = create the dummy boxes with particles in them.
#NOTE: I am not sure if addParticles is computationally demanding for other reason apart from the costs on numpy operations.
    def createDummyBox(self, location, pic, species, delta_n, n_vel, shift_vel):
        phys_loc = pic.mesh.getPosition(location)
        c = 0
        # Node by node
        for i in range(len(location)):
            # Amount of particles
            # Setting up position and velocities
            if phys_loc[i,0] == pic.mesh.xmin:
                if phys_loc[i,1] == pic.mesh.ymin:
                    mpf_new = delta_n[i]*4*pic.mesh.volumes[location[i]]/species.spwt+species.mesh_values.residuals[location[i]]+numpy.random.rand()
                    mp_new = mpf_new.astype(int)
                    species.mesh_values.residuals[location[i]] = mpf_new-mp_new
                    pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dx
                    pos[:,1] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dy
                    pos = numpy.delete(pos, numpy.flatnonzero(numpy.logical_and(pos[:,0] >= pic.mesh.xmin, pos[:,1] >= pic.mesh.ymin)), axis = 0)
                elif phys_loc[i,1] == pic.mesh.ymax:
                    mpf_new = delta_n[i]*4*pic.mesh.volumes[location[i]]/species.spwt+numpy.random.rand()
                    mp_new = mpf_new.astype(int)
                    pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dx
                    pos[:,1] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dy
                    pos = numpy.delete(pos, numpy.flatnonzero(numpy.logical_and(pos[:,0] >= pic.mesh.xmin, pos[:,1] <= pic.mesh.ymax)), axis = 0)
                else:
                    mpf_new = delta_n[i]*pic.mesh.volumes[location[i]]/species.spwt+numpy.random.rand()
                    mp_new = mpf_new.astype(int)
                    pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                    pos[:,0] -= numpy.random.rand(mp_new)*pic.mesh.dx/2
                    pos[:,1] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dy
            elif phys_loc[i,0] == pic.mesh.xmax:
                if phys_loc[i,1] == pic.mesh.ymin:
                    mpf_new = delta_n[i]*4*pic.mesh.volumes[location[i]]/species.spwt+numpy.random.rand()
                    mp_new = mpf_new.astype(int)
                    pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dx
                    pos[:,1] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dy
                    pos = numpy.delete(pos, numpy.flatnonzero(numpy.logical_and(pos[:,0] <= pic.mesh.xmax, pos[:,1] >= pic.mesh.ymin)), axis = 0)
                elif phys_loc[i,1] == pic.mesh.ymax:
                    mpf_new = delta_n[i]*4*pic.mesh.volumes[location[i]]/species.spwt+numpy.random.rand()
                    mp_new = mpf_new.astype(int)
                    pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dx
                    pos[:,1] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dy
                    pos = numpy.delete(pos, numpy.flatnonzero(numpy.logical_and(pos[:,0] <= pic.mesh.xmax, pos[:,1] <= pic.mesh.ymax)), axis = 0)
                else:
                    mpf_new = delta_n[i]*pic.mesh.volumes[location[i]]/species.spwt+numpy.random.rand()
                    mp_new = mpf_new.astype(int)
                    pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += numpy.random.rand(mp_new)*pic.mesh.dx/2
                    pos[:,1] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dy
            elif phys_loc[i,1] == pic.mesh.ymin:
                mpf_new = delta_n[i]*pic.mesh.volumes[location[i]]/species.spwt+numpy.random.rand()
                mp_new = mpf_new.astype(int)
                pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                pos[:,0] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dx
                pos[:,1] -= numpy.random.rand(mp_new)*pic.mesh.dy/2
            else:
                mpf_new = delta_n[i]*pic.mesh.volumes[location[i]]/species.spwt+numpy.random.rand()
                mp_new = mpf_new.astype(int)
                pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                pos[:,0] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dx
                pos[:,1] += numpy.random.rand(mp_new)*pic.mesh.dy/2

            vel = super().sampleIsotropicVelocity(numpy.asarray([n_vel[i]]), numpy.shape(pos)[0])+shift_vel[i,:]
            self.addParticles(species,pos,vel)

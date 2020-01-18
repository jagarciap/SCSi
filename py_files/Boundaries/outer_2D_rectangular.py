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

#	+applyElectricBoundary(Electric_Field) = Applies the boundary condition to the electric field passed as argument.
#       In this case, V = 0 at the boundaries.
    def applyElectricBoundary(self, e_field):
        pass

#	+applyMagneticBoundary(Magnetic_Field) = Applies the boundary condition to the magnetic field passed as argument.
#       No magnetic field so far
    def applyMagneticBoundary(self, m_field):
        pass

#	+applyParticleBoundary(Species) = Applies the boundary condition to the species passed as argument.
    def applyParticleBoundary(self, species):
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
        print('Number of {} eliminated:'.format(species.type), count2)

#       +createDummyBox([ind]location, PIC pic, Species species, [double] delta_n, [double] n_vel, [double] shift_vel) = create the dummy boxes with particles in them.
#NOTE: I am not sure if addParticles is computationally demanding for other reason apart from the costs on numpy operations.
    def createDummyBox(self, location, pic, species, delta_n, n_vel, shift_vel):
        phys_loc = pic.mesh.getPosition(location)
        c = 0
        # Node by node
        for i in range(len(location)):
            # Amount of particles
            mpf_new = delta_n[i]*pic.mesh.volumes[location[i]]/species.spwt+numpy.random.rand()
            mp_new = mpf_new.astype(int)
            # Setting up position and velocities
            if phys_loc[i,0] == pic.mesh.xmin:
                if phys_loc[i,1] == pic.mesh.ymin:
                    pos = numpy.ones((mp_new*4, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += (numpy.random.rand(mp_new*4)-0.5)*pic.mesh.dx
                    pos[:,1] += (numpy.random.rand(mp_new*4)-0.5)*pic.mesh.dy
                    pos = numpy.delete(pos, numpy.logical_and(pos[:,0] >= pic.mesh.xmin, pos[:,1] >= pic.mesh.ymin) , axis = 0)
                elif phys_loc[i,1] == pic.mesh.ymax:
                    pos = numpy.ones((mp_new*4, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += (numpy.random.rand(mp_new*4)-0.5)*pic.mesh.dx
                    pos[:,1] += (numpy.random.rand(mp_new*4)-0.5)*pic.mesh.dy
                    pos = numpy.delete(pos, numpy.logical_and(pos[:,0] >= pic.mesh.xmin, pos[:,1] <= pic.mesh.ymax) , axis = 0)
                else:
                    pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                    pos[:,0] -= numpy.random.rand(mp_new)*pic.mesh.dx/2
                    pos[:,1] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dy
            elif phys_loc[i,0] == pic.mesh.xmax:
                if phys_loc[i,1] == pic.mesh.ymin:
                    pos = numpy.ones((mp_new*4, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += (numpy.random.rand(mp_new*4)-0.5)*pic.mesh.dx
                    pos[:,1] += (numpy.random.rand(mp_new*4)-0.5)*pic.mesh.dy
                    pos = numpy.delete(pos, numpy.logical_and(pos[:,0] <= pic.mesh.xmax, pos[:,1] >= pic.mesh.ymin) , axis = 0)
                elif phys_loc[i,1] == pic.mesh.ymax:
                    pos = numpy.ones((mp_new*4, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += (numpy.random.rand(mp_new*4)-0.5)*pic.mesh.dx
                    pos[:,1] += (numpy.random.rand(mp_new*4)-0.5)*pic.mesh.dy
                    pos = numpy.delete(pos, numpy.logical_and(pos[:,0] <= pic.mesh.xmax, pos[:,1] <= pic.mesh.ymax) , axis = 0)
                else:
                    pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                    pos[:,0] += numpy.random.rand(mp_new)*pic.mesh.dx/2
                    pos[:,1] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dy
            elif phys_loc[i,1] == pic.mesh.ymin:
                pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                pos[:,0] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dx
                pos[:,1] -= numpy.random.rand(mp_new)*pic.mesh.dy/2
            else:
                pos = numpy.ones((mp_new, species.pos_dim))*phys_loc[i].T
                pos[:,0] += (numpy.random.rand(mp_new)-0.5)*pic.mesh.dx
                pos[:,1] += numpy.random.rand(mp_new)*pic.mesh.dy/2

            vel = super().sampleIsotropicVelocity(n_vel[i], numpy.shape(pos)[0])+shift_vel[i,:]
            self.addParticles(species,pos,vel)

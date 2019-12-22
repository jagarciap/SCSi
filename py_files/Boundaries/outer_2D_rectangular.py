import numpy

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
        e_field.potential[self.location] = 0

#	+applyMagneticBoundary(Magnetic_Field) = Applies the boundary condition to the magnetic field passed as argument.
#       No magnetic field so far
    def applyMagneticBoundary(self, m_field):
        pass

#	+applyParticleBoundary(Species) = Applies the boundary condition to the species passed as argument.
    def applyParticleBoundary(self, species):
        #Just for convenience in writing
        np = species.part_values.current_n
        #Finding the particles
        ind = numpy.flatnonzero(numpy.logical_or(numpy.logical_or(numpy.logical_or(species.part_values.position[:np,0]< self.xmin, \
                species.part_values.position[:np,0] > self.xmax), \
                species.part_values.position[:np,1] > self.ymax), \
                species.part_values.position[:np,1] < self.ymin))
        # Eliminating particles
        ndi_out = super().removeParticles(species,ind)
        count2 = numpy.shape(ind)[0]
        print('Number of{}eliminated:'.format(species.type), count2)

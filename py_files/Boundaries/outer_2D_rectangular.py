import boundary.py

#Geometrical parameters for a rectangular outer boundary
x_min = 0.0
x_max = 10.0
y_min = -5.0
y_max = 5.0


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
    def __init__(self):
        self.type = "Outer - Domain"
        self.xmin = x_min
        self.xmax = x_max
        self.ymin = y_min
        self.ymax = y_max
        location = []

#	+applyElectricBoundary(Electric_Field) = Applies the boundary condition to the electric field passed as argument.
#       In this case, V = 0 at the boundaries.
    def applyElectricBoundary(self, e_field):
        e_field.potential[location] = 0

#	+applyMagneticBoundary(Magnetic_Field) = Applies the boundary condition to the magnetic field passed as argument.
#       No magnetic field so far
    def applyMagneticBoundary(self, m_field):
        pass

#	+applyParticleBoundary(Species) = Applies the boundary condition to the species passed as argument.
    def applyParticleBoundary(self, species):
        pass

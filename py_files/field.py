#Data structures that contain the fields of the system
import numpy

import mesh as m

#Field (Abstract):
#
#Definition = Indicate the attributes and methods that all fields have to implement. The fields obtain the functions to compute the fields from 'solver.py'
#Attributes:
#	+pic (PIC) = Class that contains PIC methods.
#	+boundaries ([Boundary]) = Class that represents the methods which apply boundaries to the fields.
#	+nPoints (int) = number of nodes in the mesh (Same as mesh class).
#	+field ([double, double]) = components (columns) of field at each node.
#Methods:
#	+__init__(...) = This function, for each subclass, will take care of the initial condition of the field.
#	+computeField([Species] species) = Computes the updated field values.
#	+fieldAtParticles([double,double] position) [double,double] = return an array indicating by component (columns) and particles (rows) the field at every position.
class Field(object):
    def __init__(self, n_pic, n_boundaries, n_points, field_dim):
            self.pic = n_pic
            self.boundaries = n_boundaries
            self.nPoints = n_points
            self.field = numpy.zeros((self.nPoints, field_dim))

    def computeField(self, species):
        pass

    def fieldAtParticles(self, position):
        return self.pic.gather(position, self.field)


#Electric_Field (Inherits from Field):
#
#Definition = Electric field
#Attributes:
#	+type (string) = some string that describes the source of the electric field (created by the interaction plasma-spacecraft, user-defined, constant, etc.)
#	+potential ([double]) = Electric potential at each node of the mesh.
#	+Field attributes.
#Methods:
#	+Field methods.
class Electric_Field(Field):
    def __init__(self, n_pic, n_boundaries, n_points, field_dim, n_string):
        self.type = n_string
        self.potential = numpy.zeros((n_points))
        super().__init__(n_pic, n_boundaries, n_points, field_dim)

    def computeField(self, species):
        pass


#Constant_Electric_Field(Inherits from Electric_Field):
#
#Definition = Constant electric field impsoed by the user. Does not change through time.
#Attributes:
#	+type (string) = "Electric field - Constant".
#	+Electric_Field attributes.
#Methods:
#	+Electric_Field methods.
class Constant_Electric_Field(Electric_Field):
    def __init__(self, n_pic, n_boundaries, n_points, field_dim):
        super().__init__(n_pic, n_boundaries, n_points, field_dim, "Electric field - Constant")
        self.field[:,0] += 10.0

    def computeField(self, species):
        pass
        

#Magnetic_Field (Inherits from Field):
#
#Definition = Magnetic field
#Attibutes:
#	+type (string) = some string that describes the source of the magnetic field (created by the interaction plasma-spacecraft, user-defined, constant, etc.)
#	+Field attributes.
#Methods:
#	+Field methods.
class Magnetic_Field(Field):
    def __init__(self, n_pic, n_boundaries, n_points, field_dim, n_string):
        self.type = n_string
        super.__init__(n_pic, n_boundaries, n_points, field_dim)

    def computeField(self, species):
        pass

#Data structures that contain the fields of the system
import numpy

import constants as c
import mesh as m
import solver as slv

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
    def __init__(self, n_pic, field_dim):
            self.pic = n_pic
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
    def __init__(self, n_pic, field_dim, n_string):
        self.type = n_string
        self.potential = numpy.zeros((n_points))
        super().__init__(n_pic, field_dim)

    def computeField(self, species):
        pass


#Electrostatic_2D_rm_Electric_Field (Inherits from Electric_Field):
#
#Definition = Electric field for a 2D rectangular mesh, detached from the magnetic field. Uses methods from "solver.py" to calculate electric potential, and then electric field.
#Attributes:
#	+type (string) = "Electric field - Electrostatic_2D_rm".
#	+Elctric_Field attributes.
#Methods:
#	+Electric_Field methods.
class Electrostatic_2D_rm_Electric_Field (Electric_Field):
    def __init__(self, n_pic, field_dim):
        super().__init__(n_pic, field_dim, "Electric field - Electrostatic_2D_rm")

    def computeField(self, species):
        #Prepare the right-hand-side of the Poisson equation 
        rho = numpy.zeros_like(species[0].part_values.density)
        rho += specie.part_values.density/specie.q for specie in species
        rho /= -c.EPS_0
        slv.poissonSolver_2D_rm_SORCA(self.pic.mesh, self.potential, rho)
        self.field = -slv.derive_2D_rm(self.pic.mesh, self.potential)

#Definition = Constant electric field impsoed by the user. Does not change through time.
#Attributes:
#	+type (string) = "Electric field - Constant".
#	+Electric_Field attributes.
#Methods:
#	+Electric_Field methods.
class Constant_Electric_Field(Electric_Field):
    def __init__(self, n_pic, field_dim):
        super().__init__(n_pic, field_dim, "Electric field - Constant")
        self.field[:,0] += 0.0

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
    def __init__(self, n_pic, field_dim, n_string):
        self.type = n_string
        super.__init__(n_pic, n_boundaries, n_points, field_dim)

    def computeField(self, species):
        pass

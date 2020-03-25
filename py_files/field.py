#Data structures that contain the fields of the system
import numpy
import pdb
from vtk.util.numpy_support import vtk_to_numpy

import constants as c
import mesh as m
import solver as slv
from timing import Timing

#Field (Abstract):
#
#Definition = Indicate the attributes and methods that all fields have to implement. The fields obtain the functions to compute the fields from 'solver.py'
#Attributes:
#	+name (string) = some string that describes the source and type of the field (created by the interaction plasma-spacecraft, user-defined, constant, etc.)
#	+pic (PIC) = Class that contains PIC methods.
#	+boundaries ([Boundary]) = Class that represents the methods which apply boundaries to the fields.
#	+nPoints (int) = number of nodes in the mesh (Same as mesh class).
#	+field ([double, double]) = components (columns) of field at each node.
#Methods:
#	+__init__(...) = This function, for each subclass, will take care of the initial condition of the field.
#	+computeField([Species] species) = Computes the updated field values.
#	+fieldAtParticles([double,double] position) [double,double] = return an array indicating by component (columns) and particles (rows) the field at every position.
#       +saveVTK(Mesh mesh): dictionary = Return the attributes of field to be printed in the VTK file.
#           The process is handled inside each particular field, and the final result can be constructed from the output of different subclasses.
#       +loadVTK(Mesh mesh, output) = Takes information of the field from a VTK file through 'output' and stores it in the corresponding attributes.
#           The process is handled inside each particular field, and the final result can be constructed from the output of different subclasses.
class Field(object):
    def __init__(self, n_name, n_pic, field_dim):
        self.name = n_name
        self.pic = n_pic
        self.field = numpy.zeros((self.pic.mesh.nPoints, field_dim))

    def computeField(self, species):
        pass

    def fieldAtParticles(self, position):
        return self.pic.gather(position, self.field)

    def saveVTK(self, mesh):
        return {self.name+"-field" : mesh.vtkOrdering(self.field)}

    def loadVTK(self, mesh, output):
        self.field = mesh.reverseVTKOrdering(vtk_to_numpy(output.GetPointData().GetArray(self.name+"-field")))


#Electric_Field (Inherits from Field):
#
#Definition = Electric field
#Attributes:
#	+potential ([double]) = Electric potential at each node of the mesh.
#	+Field attributes.
#Methods:
#	+Field methods.
class Electric_Field(Field):
    def __init__(self, n_pic, field_dim, n_string):
        self.potential = numpy.zeros((n_pic.mesh.nPoints))
        super().__init__("Electric"+n_string, n_pic, field_dim)

    def computeField(self, species):
        pass

    def saveVTK(self, mesh):
        dic = super().saveVTK(mesh)
        dic[self.name+"-potential"] = mesh.vtkOrdering(self.potential)
        return dic

    def loadVTK(self, mesh, output):
        self.potential = mesh.reverseVTKOrdering(vtk_to_numpy(output.GetPointData().GetArray(self.name+"-potential")))
        super().loadVTK(mesh, output)


#Electrostatic_2D_rm_Electric_Field (Inherits from Electric_Field):
#
#Definition = Electric field for a 2D rectangular mesh, detached from the magnetic field. Uses methods from "solver.py" to calculate electric potential, and then electric field.
#Attributes:
#	+type (string) = "Electric - Electrostatic_2D_rm".
#	+Elctric_Field attributes.
#Methods:
#	+Electric_Field methods.
class Electrostatic_2D_rm(Electric_Field):
    def __init__(self, n_pic, field_dim):
        super().__init__(n_pic, field_dim, " - Electrostatic_2D_rm")

    @Timing
    def computeField(self, species):
        #Prepare the right-hand-side of the Poisson equation 
        rho = numpy.zeros_like(species[0].mesh_values.density)
        for specie in species:
            rho += specie.mesh_values.density*specie.q
        rho /= -c.EPS_0
        #pdb.set_trace()
        slv.poissonSolver_2D_rm_SORCA(self.pic.mesh, self.potential, rho)
        self.field = -slv.derive_2D_rm(self.pic.mesh, self.potential)
        for boundary in self.pic.mesh.boundaries:
            boundary.applyElectricBoundary(self)

#       +Computation of Dirichlet boundary condition at every node in location ([ind]). Every row in value ([double]) corresponds to one node in location.
    def dirichlet(self, location, values):
        #Dirichlet
        self.potential[location] = values
        #Electric field trough Pade 2nd order in the boundaries
        self.field[location, :] = -slv.derive_2D_rm_boundaries(location, self.pic.mesh, self.potential)

#       +Neumann([ind] location, [double] valus) = Set Neumann conditions in the nodes at 'location'.
#           +values account for the values of the e_field normal to the border.
#           +Note: The Function doesn't handle the situation of the corners.
    def neumann(self, location, values):
        # Variables at hand
        nx = self.pic.mesh.nx
        ny = self.pic.mesh.ny
        dx = self.pic.mesh.dx
        dy = self.pic.mesh.dy
        #Field and potential
        for i in range(len(location)):
            if location[i] < nx:
                self.field[location[i],1] = values[i]
                self.potential[location[i]] = self.field[location[i],1]*dy+self.potential[location[i]+nx]
            elif location[i] > nx*(ny-1):
                self.field[location[i],1] = values[i]
                self.potential[location[i]] = -self.field[location[i],1]*dy+self.potential[location[i]-nx]
            elif location[i]%nx == 0:
                self.field[location[i],0] = values[i]
                self.potential[location[i]] = self.field[location[i],0]*dx+self.potential[location[i]+1]
            else:
                self.field[location[i],0] = values[i]
                self.potential[location[i]] = -self.field[location[i], 0]*dx+self.potential[location[i]-1]
            
#Constant_Electric_Field(Electric_Field):
#
#Definition = Constant electric field impsoed by the user. Does not change through time.
#Attributes:
#	+type (string) = "Electric field - Constant".
#	+Electric_Field attributes.
#Methods:
#	+Electric_Field methods.
class Constant_Electric_Field(Electric_Field):
    def __init__(self, n_pic, field_dim):
        super().__init__(n_pic, field_dim, " - Constant")
        self.field[:,0] += 0.27992

    def computeField(self, species):
        pass
        
#Time_Electric_Field(Electric_Field):
#
#Definition = Electric field dependent on time.
#Attributes:
#	+type (string) = "Electric field - Constant".
#	+Electric_Field attributes.
#Methods:
#	+Electric_Field methods.
#       +computeField(Species species, int p_step, int e_step = 0) = Recieves the steps in the simulation and computes the time from the start of the execution. Then, updates the field accordingly
#           with this and any function imposed by the user inside of this method.
class Time_Electric_Field(Electric_Field):
    def __init__(self, n_pic, field_dim):
        super().__init__(n_pic, field_dim, " - Constant")
        #Initial electric field
        self.field[:,0] += 0.27992

    def computeField(self, species, p_step, e_step = 0):
        time = p_step*c.P_DT+e_step*c.E_DT
        gradient = 98.4658
        self.field[:,0] = 0.27992+gradient*time

#Magnetic_Field (Inherits from Field):
#
#Definition = Magnetic field
#Attibutes:
#	+Field attributes.
#Methods:
#	+Field methods.
class Magnetic_Field(Field):
    def __init__(self, n_pic, field_dim, n_string):
        super().__init__("Magnetic"+n_string, n_pic, field_dim)

    def computeField(self, species):
        pass

#Constant_Magnetic_Field(Inherits from Magnetic_Field):
#
#Definition = Constant Magnetic field impsoed by the user. Does not change through time. It works as a perpendicular field to the 2D dominion of the electric field and particular. Thus, field_dim = 1.
#Attributes:
#	+type (string) = "Electric - Constant".
#	+Magnetic_Field attributes.
#Methods:
#	+Magnetic_Field methods.
class Constant_Magnetic_Field(Magnetic_Field):
    def __init__(self, n_pic, field_dim):
        super().__init__(n_pic, field_dim, " - Constant")
        self.field += c.B_STRENGTH

    def computeField(self, species):
        pass

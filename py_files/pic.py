#Data structures that contain PIC numerical methods
import numpy

import constants as c
import mesh as m

#PIC (Abstract)(Association between Mesh and PIC):
#
#Definition = Indicate the methods that all PIC classes have to implement. Each PIC concrete class will depend on the type of mesh, as well as the type of PIC algorithm implemented.
#Attributes:
#	+mesh (Mesh) = Instance of the mesh class for later use of getIndex().
#Methods:
#	+scatter([double, double] positions, [double] values, [double] field) = Receives the positions of the particles, and makes scatter procedure, calculating the values of field for each node.
#	+gather([double, double] positions, [double, double] field): [double, double]field_p = Calculates values of the field in particles' positions, returning these values in an array as long as positions,
#           The columns are the (x,y,z) components
#       +scatterDiffSq([double, double] positions, [double] values, [double] array_diff, [double] field) = Makes a PIC averaging over 
#           (values-array_diff)**2 in all the nodes involved for every particle in positions. 
#           Thus, array_diff an field need to have the same dimension, and values have to be broadcastable to the len of positions.
class PIC (object):
    def __init__(self, mesh):
        self.mesh = mesh

    def scatter(self, positions, values, field):
        pass

    def scatterDiffSq(self, positions, values, array_diff, field):
        pass

    def gather(self, positions, field):
        pass


#PIC_2D_rm1o (Inherits from PIC):
#
#Definition = PIC class for rm10 ('r'ectangular 'm'esh, '1'st 'o'rder implementation).
#Attributes:
#	+PIC class attributes.
#Methods:
#	+Implementation of PIC class methods.
#	+scatterDensity (Species) = return densities of that species in every node of the mesh.
#	+scatterVelocity (Species) = return velocities of that species in every node of the mesh.
#	+scatterFlux = return flux of particles of that species into every indicated node (not all the mesh).
class PIC_2D_rm1o(PIC):

    string = "PIC method for a rectangular 2D mesh, making interpolations at first order"

    def __init__(self, mesh):
        super(PIC_2D_rm1o, self).__init__(mesh)

    def scatter(self, positions, values, field):
        #Getting mesh coordinates
        mc = self.mesh.getIndex(positions)
        
        index = mc.astype(int)
        array = self.mesh.indexToArray(index)
        di = mc[:,0] - index[:,0]
        dj = mc[:,1] - index[:,1]
        
        numpy.add.at(field, array, (1-di)*(1-dj)*values)
        numpy.add.at(field, array+1, di*(1-dj)*values)
        numpy.add.at(field, array+self.mesh.nx, (1-di)*dj*values)
        numpy.add.at(field, array+self.mesh.nx+1, di*dj*values)

#       +scatterDiffSq([double, double] positions, [double] values, [double] array_diff, [double] field) = Makes a PIC averaging over 
#           (values-array_diff)**2 in all the nodes involved for every particle in positions. 
#           Thus, array_diff an field need to have the same dimension, and values have to be broadcastable to the len of positions.
    def scatterDiffSq(self, positions, values, array_diff, field):
        #Getting mesh coordinates
        mc = self.mesh.getIndex(positions)
        
        index = mc.astype(int)
        array = self.mesh.indexToArray(index)
        di = mc[:,0] - index[:,0]
        dj = mc[:,1] - index[:,1]
        
        numpy.add.at(field, array, (1-di)*(1-dj)*(values-array_diff[array])*(values-array_diff[array]))
        numpy.add.at(field, array+1, di*(1-dj)*(values-array_diff[array+1])*(values-array_diff[array+1]))
        numpy.add.at(field, array+self.mesh.nx, (1-di)*dj*(values-array_diff[array+self.mesh.nx])*(values-array_diff[array+self.mesh.nx]))
        numpy.add.at(field, array+self.mesh.nx+1, di*dj*(values-array_diff[array+self.mesh.nx+1])*(values-array_diff[array+self.mesh.nx+1]))

#       +scatterDensity (Species) = return densities of that species in every node of the mesh.
    def scatterDensity(self, species):
        #reset the density
        species.mesh_values.density *= 0
    
        #scatter particles to the mesh
        self.scatter(species.part_values.position[:species.part_values.current_n], species.spwt, species.mesh_values.density)
        
        #divide by cell volume
        species.mesh_values.density /= self.mesh.volumes
    
#       +scatterVelocity (Species) = return velocities of that species in every node of the mesh.
    def scatterSpeed(self, species):
        #reset the velocity
        species.mesh_values.velocity *= 0
    
        #scatter particles to the mesh
        for dim in range(numpy.shape(species.part_values.velocity)[1]):
            self.scatter(species.part_values.position[:species.part_values.current_n], \
                    species.part_values.velocity[:species.part_values.current_n,dim], species.mesh_values.velocity[:,dim])
            #Finalizing velocity
            species.mesh_values.velocity[:,dim] *= numpy.where(species.mesh_values.density < 1e-5, 0.0, species.spwt/species.mesh_values.density/self.mesh.volumes)

#       +scatterTemperature(Species) = return temperature of that species in every node of the mesh.
    def scatterTemperature(self, species):
        #Reset temperature
        species.mesh_values.temperature *= 0
        #Scatter temperature
        for dim in range(numpy.shape(species.part_values.velocity)[1]):
            self.scatterDiffSq(species.part_values.position[:species.part_values.current_n], \
                    species.part_values.velocity[:species.part_values.current_n,dim], species.mesh_values.velocity[:,dim], species.mesh_values.temperature)
        #Finalizing temperature
        species.mesh_values.temperature *= numpy.where(species.mesh_values.density < 1e-5, 0.0, species.spwt/species.mesh_values.density/self.mesh.volumes*species.m/c.K)
    
#       +scatterFlux = return flux of particles of that species into every indicated node (not all the mesh).
    def scatterFlux(self):
        pass
    
#	+gather([double, double] positions, [double, double] field): [double, double]field_p = Calculates values of the field in particles' positions, returning these values in an array as long as positions,
#                                                                                               The columns are the (x,y,z) positions
    def gather(self, positions, field):
        #Getting mesh coordinates
        mc = self.mesh.getIndex(positions)
        #Creating the array
        values = numpy.zeros((numpy.shape(positions)[0], numpy.shape(field)[1]))
        
        index = mc.astype(int)
        array = self.mesh.indexToArray(index)
        di = mc[:,0] - index[:,0]
        dj = mc[:,1] - index[:,1]

        #NOTE: Maybe this can be further optmized later
        di = numpy.repeat(di[:,None], 2, axis = 1)
        dj = numpy.repeat(dj[:,None], 2, axis = 1)

        #From mesh to particles, summing in every dimension through different columns
        values += field[array,:]*(1-di)*(1-dj)
        values += field[array+1,:]*di*(1-dj)
        values += field[array+self.mesh.nx,:]*(1-di)*dj
        values += field[array+self.mesh.nx+1,:]*di*dj

        return values

import numpy

#Motion_Solver(Abstract-like: As with Boundary, some methods can be applied to various motion solvers, so they will be stored here):
#
#Definition = Class that shows which methods and attributes need to have all motion solvers.
#Attributes:
#	+type (string) = string indicating the type of scheme.
#       +pic (PIC) = PIC solver.
#Methods:
#       +initialConfiguration(Species, [Field]) = Make necessary adjustment to the initial configuration of species so as to begin the advancement in time.
#	+advance(Species, [Field]) = Advance the particles in time. It will treat the particles as well as update the mesh_values.
#	+updateMeshValues(Species) = Update the attributes of Particles_In_Mesh.
#       +updateParticles(Species, [Field]) = Particle advance in time.
class Motion_Solver(object):
    def __init__(self, pic_slv):
        self.pic = pic_slv

    def initialConfiguration(self, species, fields):
        pass

    def advance(self, species, fields):
        pass
    
    def updateMeshValues(self, species):
        pass

    def updateParticles(self, species, fields):
        pass


#Leap_Frog(Inherits from Motion_Solver):
#
#Definition = Implementation of Leap_Frog method. This method, as being more specific, is more dependant on the current situation of each simulation. As such it has to
#               be updated or new classes has to be created if the situation varies.
#
#Attributes: 
#	+type (string) = "Leap Frog"
#       +pic (PIC) = PIC solver. For this class specific methods not provided by PIC super class are used.
#Methods:
#       +initialConfiguration(Species, Field) = Make necessary adjustment to the initial configuration of species so as to begin the advancement in time.
#           it also takes care of the first update of Particle values in mesh. So far just E, so [Field]->Field.
#       +rewindVelocity(species, field) = Take the velocity of particles half a step back in time for 'field' electric field.
#       +updateParticles(Species, Field) = Particle advance in time. So far only E, so [Field]->Field in argument.
#	+Motion_Solver methods.
class Leap_Frog(Motion_Solver):
    def __init__(self, pic_slv):
        super().__init__(pic_slv)

#       +initialConfiguration(Species, Field) = Make necessary adjustment to the initial configuration of species so as to begin the advancement in time.
#           So far just E, so [Field]->Field.
    def initialConfiguration(self, species, field):
        self.rewindVelocity(species, field)

#       +rewindVelocity(species, field) = Take the velocity of particles half a step back in time for 'field' electric field.
    def rewindVelocity(self, species, field):
        species.part_values.velocity[:species.part_values.current_n,:] -= species.q_over_m*species.dt/2*field.fieldAtParticles(species.part_values.position[:species.part_values.current_n])

#	+advance(Species, [Field]) = Advance the particles in time. It will treat the particles as well as update the mesh_values.
    def advance(self, species, fields):
        self.updateParticles(species, fields)
        self.updateMeshValues(species)

#	+updateMeshValues(Species) = Update the attributes of Particles_In_Mesh. Particular for each species, so it needs to be updated with every new species.
    def updateMeshValues(self, species):
        if species.name == "Electron - Solar wind" or species.name == "Proton - Solar wind":
            self.pic.scatterDensity(species)
            self.pic.scatterSpeed(species)

#       +updateParticles(Species, Field) = Particle advance in time. So far only E, so [Field]->Field in argument.
    def updateParticles(self, species, field):
        np = species.part_values.current_n
        species.part_values.velocity[:np,:] += species.q_over_m*species.dt*field.fieldAtParticles(species.part_values.position[:np,:])
        species.part_values.position[:np,:] += species.part_values.velocity[:np,:]*species.dt
        self.pic.mesh.boundaries[0].applyParticleBoundary(species)

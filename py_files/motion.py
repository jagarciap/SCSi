import copy
from functools import reduce
import numpy
import pdb

from field import Field
from timing import Timing

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
#       +vel_dic {String : [double,double]} = Is a dictionary containing the new values of velocities for every species. The key 'String' is the actual name of the species.
#Methods:
#       +initialConfiguration(Species, Field) = Make necessary adjustment to the initial configuration of species so as to begin the advancement in time.
#           So far just E, so [Field]->Field.
#       +rewindVelocity(species, field) = Take the velocity of particles half a step back in time for 'field' electric field.
#	+advance(Species, [Field]) = Advance the particles in time. It will treat the particles as well as update the mesh_values.
#           Extent is a karg for updateMeshValues().
#	+updateMeshValues(Species) = Update the attributes of Particles_In_Mesh. Particular for each species, so it needs to be updated with every new species.
#           Extent can be '0' indicating that every attribute of mesh_values is updated, '1' indicating that only necessary attributes are updated, and 
#           '2' indicating that the 'non-necessary' attributes are the only ones updated. Notice that the criteria between 'necessary' and 'non-necessary'
#           attributes will change with the type of phenomena being included.
#       +updateParticles(Species, Field) = Particle advance in time. So far only E, so [Field]->Field in argument.
#       +motionTreatment(Species species) = It takes the velocities array in the species.part_values atribute and average it with the velocities stored in vel_dic for the particular species.
#           That is, to synchronize velocity with position.
#	+Motion_Solver methods.
class Leap_Frog(Motion_Solver):
    def __init__(self, pic_slv, species_names, max_n, vel_dim):
        super().__init__(pic_slv)
        self.type = "Leap Frog"
        self.vel_dic = {}
        for i in range(len(species_names)):
            self.vel_dic[species_names[i]] = numpy.array((max_n[i], vel_dim[i]))

#       +initialConfiguration(Species, Field) = Make necessary adjustment to the initial configuration of species so as to begin the advancement in time.
#           So far just E, so [Field]->Field.
    def initialConfiguration(self, species, field):
        self.rewindVelocity(species, field)

#       +rewindVelocity(species, field) = Take the velocity of particles half a step back in time for 'field' electric field.
    def rewindVelocity(self, species, field):
        species.part_values.velocity[:species.part_values.current_n,:] -= species.q_over_m*species.dt/2*field.fieldAtParticles(species.part_values.position[:species.part_values.current_n])

#	+advance(Species, [Field]) = Advance the particles in time. It will treat the particles as well as update the mesh_values.
#           extent is a karg for updateMeshValues().
#           update_dic is a karg for updateParticles(). 
#           type_boundary indicates the type of boundary method to apply to particles. 'open', the default mthod, deletes them. 'reflective' reflects them back to the dominion.
#           **kwargs may contain arguments necessary for inner methods.
    @Timing
    def advance(self, species, fields, extent = 0, update_dic = 1, type_boundary = 'open', **kwargs):
        self.updateParticles(species, fields, update_dic, type_boundary, **kwargs)
        self.updateMeshValues(species, extent)

#	+updateMeshValues(Species) = Update the attributes of Particles_In_Mesh. Particular for each species, so it needs to be updated with every new species.
#           Extent can be '0' indicating that every attribute of mesh_values is updated, '1' indicating that only necessary attributes are updated, and 
#           '2' indicating that the 'non-necessary' attributes are the only ones updated. Notice that the criteria between 'necessary' and 'non-necessary'
#           attributes will change with the type of phenomena being included.
    @Timing
    def updateMeshValues(self, species, extent = 0):
        #if species.name == "Electron - Solar wind" or species.name == "Proton - Solar wind":
        if extent == 0:
            self.pic.scatterDensity(species)
            self.motionTreatment(species)
            self.pic.scatterSpeed(species)
            self.pic.scatterTemperature(species)
        elif extent == 1:
            self.pic.scatterDensity(species)
        elif extent == 2:
            self.motionTreatment(species)
            self.pic.scatterSpeed(species)
            self.pic.scatterTemperature(species)

#       +updateParticles(Species, Field, int) = Particle advance in time. So far only E, so [Field]->Field in argument.
#           If update_dic == 1, the vel_dic entry for the species is updated, otherwise it remains untouched.
#           type_boundary indicates the type of boundary method to apply to particles. 'open', the default mthod, deletes them. 'reflective' reflects them back to the dominion.
#           **kwargs may contain arguments necessary for inner methods.
    @Timing
    def updateParticles(self, species, field, update_dic, type_boundary, **kwargs):
        np = species.part_values.current_n
        species.part_values.velocity[:np,:] += species.q_over_m*species.dt*field.fieldAtParticles(species.part_values.position[:np,:])
        species.part_values.position[:np,:] += species.part_values.velocity[:np,:]*species.dt
        if update_dic == 1:
            self.vel_dic[species.name] = copy.copy(species.part_values.velocity)
        for boundary in self.pic.mesh.boundaries:
            boundary.applyParticleBoundary(species, type_boundary, **kwargs)

#       +motionTreatment(Species species) = It takes the velocities array in the species.part_values atribute and average it with the velocities stored in vel_dic for the particular species.
#           That is, to synchronize velocity with position.
    def motionTreatment(self, species):
        species.part_values.velocity[:species.part_values.current_n,:] = (species.part_values.velocity[:species.part_values.current_n,:]+\
                self.vel_dic[species.name][:species.part_values.current_n,:])/2

#Boris_Push(Inherits from Motion_Solver):
#
#Definition = Implementation of  Boris algorithm for time integration. This method also creates desynchronization between v and x.
#
#Attributes: 
#	+type (string) = "Boris Push"
#       +pic (PIC) = PIC solver. For this class specific methods not provided by PIC super class are used.
#       +vel_dic {String : [double,double]} = Is a dictionary containing the new values of velocities for every species. The key 'String' is the actual name of the species.
#Methods:
#       +initialConfiguration(Species, Field) = Make necessary adjustment to the initial configuration of species so as to begin the advancement in time.
#           So far just E, so [Field]->Field.
#       +rewindVelocity(species, field) = Take the velocity of particles half a step back in time for 'field' electric field.
#	+advance(Species, [Field]) = Advance the particles in time. It will treat the particles as well as update the mesh_values.
#           Extent is a karg for updateMeshValues().
#	+updateMeshValues(Species) = Update the attributes of Particles_In_Mesh. Particular for each species, so it needs to be updated with every new species.
#           Extent can be '0' indicating that every attribute of mesh_values is updated, '1' indicating that only necessary attributes are updated, and 
#           '2' indicating that the 'non-necessary' attributes are the only ones updated. Notice that the criteria between 'necessary' and 'non-necessary'
#           attributes will change with the type of phenomena being included.
#       +updateParticles(Species, Field) = Particle advance in time. So far only E, so [Field]->Field in argument.
#       +motionTreatment(Species species) = It takes the velocities array in the species.part_values atribute and average it with the velocities stored in vel_dic for the particular species.
#           That is, to synchronize velocity with position.
#	+Motion_Solver methods.
class Boris_Push(Motion_Solver):
    def __init__(self, pic_slv, species_names, max_n, vel_dim):
        super().__init__(pic_slv)
        self.type = "Leap Frog - Boris Push"
        self.vel_dic = {}
        for i in range(len(species_names)):
            self.vel_dic[species_names[i]] = numpy.array((max_n[i], vel_dim[i]))

#       +initialConfiguration(Species, Field) = Make necessary adjustment to the initial configuration of species so as to begin the advancement in time.
#           So far just E, so [Field]->Field. As far as I know there is no need to make a 'rewind' due to magnetic field.
    def initialConfiguration(self, species, field):
        self.rewindVelocity(species, field)

#       +rewindVelocity(species, field) = Take the velocity of particles half a step back in time for 'field' electric field.
    def rewindVelocity(self, species, field):
        species.part_values.velocity[:species.part_values.current_n,:] -= species.q_over_m*species.dt/2*field.fieldAtParticles(species.part_values.position[:species.part_values.current_n])

#	+advance(Species, [Field]) = Advance the particles in time. It will treat the particles as well as update the mesh_values.
#           extent is a karg for updateMeshValues().
#           update_dic is a karg for updateParticles(). 
#           type_boundary indicates the type of boundary method to apply to particles. 'open', the default mthod, deletes them. 'reflective' reflects them back to the dominion.
#           **kwargs may contain arguments necessary for inner methods.
    @Timing
    def advance(self, species, e_fields, m_fields, extent = 0, update_dic = 1, type_boundary = 'open', **kwargs):
        self.updateParticles(species, e_fields, m_fields, update_dic, type_boundary, **kwargs)
        self.updateMeshValues(species, extent)

#	+updateMeshValues(Species) = Update the attributes of Particles_In_Mesh. Particular for each species, so it needs to be updated with every new species.
#           Extent can be '0' indicating that every attribute of mesh_values is updated, '1' indicating that only necessary attributes are updated, and 
#           '2' indicating that the 'non-necessary' attributes are the only ones updated. Notice that the criteria between 'necessary' and 'non-necessary'
#           attributes will change with the type of phenomena being included.
    @Timing
    def updateMeshValues(self, species, extent = 0):
        #if species.name == "Electron - Solar wind" or species.name == "Proton - Solar wind":
        if extent == 0:
            self.pic.scatterDensity(species)
            self.motionTreatment(species)
            self.pic.scatterSpeed(species)
            self.pic.scatterTemperature(species)
        elif extent == 1:
            self.pic.scatterDensity(species)
        elif extent == 2:
            self.motionTreatment(species)
            self.pic.scatterSpeed(species)
            self.pic.scatterTemperature(species)

#       +updateParticles(Species, [Field], int) = Particle advance in time. So far only one E and B.
#           If update_dic == 1, the vel_dic entry for the species is updated, otherwise it remains untouched.
#           type_boundary indicates the type of boundary method to apply to particles. 'open', the default method, deletes them. 'reflective' reflects them back to the dominion.
#           **kwargs may contain arguments necessary for inner methods.
    @Timing
    def updateParticles(self, species, e_fields, m_fields, update_dic, type_boundary, **kwargs):
        #Summing over the different fields
        e_total = copy.deepcopy(e_fields[0])
        m_total = copy.deepcopy(m_fields[0])
        e_total.field = reduce(lambda x,y : x.field + y.field, e_fields)
        m_total.field = reduce(lambda x,y : x.field + y.field, m_fields)
        #Updating velocity
        np = species.part_values.current_n
        species.part_values.velocity[:np,:] = self.electricAdvance(e_total.fieldAtParticles(species.part_values.position[:np,:]), species, species.dt/2)
        species.part_values.velocity[:np,:] = self.magneticRotation(b_total.fieldAtParticles(species.part_values.position[:np,:]), species, species.dt)
        species.part_values.velocity[:np,:] = self.electricAdvance(e_total.fieldAtParticles(species.part_values.position[:np,:]), species, species.dt/2)
        #Updating position         
        species.part_values.position[:np,:] += species.part_values.velocity[:np,:]*species.dt
        if update_dic == 1:
            self.vel_dic[species.name] = copy.copy(species.part_values.velocity)
        for boundary in self.pic.mesh.boundaries:
            boundary.applyParticleBoundary(species, type_boundary, **kwargs)

#       +electricAdvance(Field field, Species species, double dt) = Updates the velocity one 'dt' forward.
#           Boris (1970) advancement in time. Is the usual Leap-Frog advancement in time.
    def electricAdvance(self, e_field, species, dt):
        return species.part_values.velocity[:species.part_values.current_n,:]+e_field*species.q_over_m*dt
    
#       +magneticRotation (Field B, Species species, double dt) = Update velocity applying magnetic rotation.
#           Buneman (1973) rotation.
    def magneticRotation (self, B, species, dt):
        new = numpy.zeros_like(species.part_values.velocity)
        t = species.q_over_m*B*dt/2
        sine = 2*t/(1+t*t)
        v_1 = species.part_values.velocity[:,0] + species.part_values.velocity[:,1]*t
        new[:,1] = species.part_values.velocity[:,1] - v_1*sine
        new[:,0] = v_1 + new[:,1]*t
        return new

#       +motionTreatment(Species species) = It takes the velocities array in the species.part_values atribute and average it with the velocities stored in vel_dic for the particular species.
#           That is, to synchronize velocity with position.
    def motionTreatment(self, species):
        species.part_values.velocity[:species.part_values.current_n,:] = (species.part_values.velocity[:species.part_values.current_n,:]+\
                self.vel_dic[species.name][:species.part_values.current_n,:])/2

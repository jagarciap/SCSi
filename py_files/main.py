# Main file of the program
import copy
import numpy
import pdb
import sys
numpy.seterr(invalid='ignore', divide='ignore', over = 'raise')

import constants as c
from field import Electrostatic_2D_rm, Constant_Electric_Field
from mesh import Mesh_2D_rm
from Species.proton import Proton_SW
from Species.electron import Electron_SW
from Species.user_defined import User_Defined
from pic import PIC_2D_rm1o
from motion import Leap_Frog
import output as out



## ---------------------------------------------------------------------------------------------------------------
# Initiating data structures for the system
## ---------------------------------------------------------------------------------------------------------------


#System:
#
#Definition = Is the class that contains every variable and class necessary for the simulation to be executed.
#Attributes:
#	+ts (int) = Timestep of the simulation.
#	+The rest of the variables will change with the simulation, but normally, there are:
#	+mesh (Mesh).
#	+pic (PIC).
#	+fields (Field) = Probably several of them.
#	+species (species) = Probably several of them.
#	+part_solver (Motion_Solver).
#Methods:
#	+Remark about init(): It will "declare" the attributes necessary for the simulation to run. The actual assignment of atributes
#		to instances of each class will occur during the 'initial condition' section of 'main.py'.
#	+arrangePickle() : Variable = Return a tuple of keys to iterate over when saving and loading of/from '.pkl' files, in the order required.
#	+arrangeVTK() : Variable = Return a tuple of keys to iterate over when saving and loading of/from VTK files, in the order required.
class System(object):
    def __init__(self):
        self.at = {}
        self.at['ts'] = 0
        self.at['mesh'] = Mesh_2D_rm(c.XMIN, c.XMAX, c.YMIN, c.YMAX, c.NX, c.NY, c.DEPTH)
        self.at['pic'] = PIC_2D_rm1o(self.at['mesh'])
        self.at['electrons'] = Electron_SW(0.0, c.E_SPWT, c.E_SIZE, c.DIM, c.DIM, self.at['mesh'].nPoints, c.NUM_TRACKED)
        self.at['protons'] = Proton_SW(0.0, c.P_SPWT, c.P_SIZE, c.DIM, c.DIM, self.at['mesh'].nPoints, c.NUM_TRACKED)
        self.at['user'] = User_Defined(c.P_DT, -c.QE, c.MP, 0, c.P_SPWT, 1, c.DIM, c.DIM, self.at['mesh'].nPoints, 0, "1")
        self.at['e_field'] = Electrostatic_2D_rm(self.at['pic'], c.DIM)
        self.at['part_solver'] = Leap_Frog(self.at['pic'], [self.at['electrons'].name, self.at['protons'].name],\
                [self.at['electrons'].part_values.max_n, self.at['protons'].part_values.max_n],\
                [self.at['electrons'].vel_dim, self.at['protons'].vel_dim])

    def arrangePickle(self):
        return ('ts', 'e_field', 'electrons', 'protons', 'user', 'part_solver')

    def arrangeVTK(self):
        return ('ts', 'e_field', 'electrons', 'protons', 'user')

    def arrangeParticlesTXT(self):
        return ('ts', 'electrons', 'protons', 'user')

#Initialization of the system
system = System()

## ---------------------------------------------------------------------------------------------------------------
# Initial condition
## ---------------------------------------------------------------------------------------------------------------
# Initial conditions are selected with a number next to the file when the program is executed, e.g. 'python3 main.py 2'.
# If the initial condition requires more arguments, they are written here in the file.
# Listed Initial conditions:
# 1: Basic empty system.
# 2: Execution from a VTK file.
# 3: Execution from a Pickle file.

if sys.argv[1] == '1':
    system.at['part_solver'].initialConfiguration(system.at['electrons'], system.at['Field'])
    system.at['part_solver'].initialConfiguration(system.at['protons'], system.at['Field'])

elif sys.argv[1] == '2':
    #File to be used as source of initial condition
    filename = 'ts00090.vtr'
    out.loadVTK(filename, system.at['mesh'], system.at, system.arrangeVTK())
    system.at['part_solver'].initialConfiguration(system.at['electrons'], system.at['Field'])
    system.at['part_solver'].initialConfiguration(system.at['protons'], system.at['Field'])

elif sys.argv[1] == '3':
    #File to be used as source of initial condition
    filename = 'sys_ts=2000_2020-03-06_21h05m.pkl'
    out.loadPickle(filename, system.at, system.arrangePickle())
    system.at['part_solver'] = Leap_Frog(system.at['pic'], [system.at['electrons'].name, system.at['protons'].name],\
            [system.at['electrons'].part_values.max_n, system.at['protons'].part_values.max_n],\
            [system.at['electrons'].vel_dim, system.at['protons'].vel_dim])

else:
    raise("Somehing is wrong here")

#Initialization of the previous step
old_system = copy.deepcopy(system)

## ---------------------------------------------------------------------------------------------------------------
# Set up of the system before the Main loop
## ---------------------------------------------------------------------------------------------------------------

#Injection of particles at all outer boundaries
#Electrons
e_n = c.E_N*numpy.ones((len(system.at['pic'].mesh.boundaries[0].location)))
drift_e_vel = numpy.zeros((len(system.at['pic'].mesh.boundaries[0].location),c.DIM))
thermal_e_vel = c.E_V_TH*numpy.ones((len(system.at['pic'].mesh.boundaries[0].location)))
#drift_e_vel[:,0] += c.E_V_SW
#Protons
p_n = c.P_N*numpy.ones((len(system.at['pic'].mesh.boundaries[0].location)))
drift_p_vel = numpy.zeros((len(system.at['pic'].mesh.boundaries[0].location),c.DIM))
thermal_p_vel = c.P_V_TH*numpy.ones((len(system.at['pic'].mesh.boundaries[0].location)))
#drift_p_vel[:,0] += c.P_V_SW
#User defined
system.at['user'] = User_Defined(c.P_DT, -c.QE, c.MP, 0, c.P_SPWT, 100, c.DIM, c.DIM, system.at['mesh'].nPoints, 0, "1")
vel = numpy.zeros((100,2))
pos = numpy.zeros((100,2))
pos[:,0] = 1.0
system.at['mesh'].boundaries[0].addParticles(system.at['user'], pos, vel)
system.at['part_solver'].updateMeshValues(system.at['user'], extent = 1)

for boundary in system.at['mesh'].boundaries:
    boundary.injectParticlesDummyBox(boundary.location, system.at['part_solver'], system.at['e_field'], system.at['electrons'], e_n, thermal_e_vel, drift_e_vel)
    boundary.injectParticlesDummyBox(boundary.location, system.at['part_solver'], system.at['e_field'], system.at['protons'], p_n, thermal_p_vel, drift_p_vel)

#Update of mesh values
system.at['part_solver'].updateMeshValues(system.at['protons'], extent = 1)
system.at['part_solver'].updateMeshValues(system.at['electrons'], extent = 1)

## ---------------------------------------------------------------------------------------------------------------
# Main loop
## ---------------------------------------------------------------------------------------------------------------

# try/except block to capture any error and save before the step before the crash. It also allows the simulation to be stopped.
try:
    while system.at['ts'] < c.NUM_TS:
        print('ts = ', system.at['ts'])
    
        # Electron motion
        for te in range(c.ELECTRON_TS):
            print('te = ', te)
            system.at['e_field'].computeField([system.at['protons'], system.at['electrons'], system.at['user']])
            system.at['part_solver'].advance(system.at['electrons'], system.at['e_field'], extent = 1)
            # Applying field borders and injecting at['electrons']
            for boundary in system.at['mesh'].boundaries:
                boundary.injectParticlesDummyBox(boundary.location, system.at['part_solver'], system.at['e_field'], system.at['electrons'], e_n, thermal_e_vel, drift_e_vel)
    
        #Proton motion
        system.at['part_solver'].advance(system.at['protons'], system.at['e_field'], extent = 1)
        for boundary in system.at['mesh'].boundaries:
            boundary.injectParticlesDummyBox(boundary.location, system.at['part_solver'], system.at['e_field'], system.at['protons'], p_n, thermal_p_vel, drift_p_vel)
    
        #Output vtk
        if system.at['ts']%10 == 0:
            system.at['part_solver'].updateMeshValues(old_system.at['electrons'], extent = 2)
            system.at['part_solver'].updateMeshValues(old_system.at['protons'], extent = 2)
            out.saveVTK(system.at['mesh'], old_system.at, system.arrangeVTK())
        if system.at['ts']%20 == 0:
            out.saveParticlesTXT(old_system.at, system.arrangeParticlesTXT())
        if system.at['ts']%1 == 0:
            out.particleTracker(old_system.at['ts'], old_system.at['protons'], old_system.at['electrons'])
    
        #Updating previous state
        old_system = copy.deepcopy(system)
    
        #Advance in timestep
        system.at['ts'] += 1

except KeyboardInterrupt:
    out.savePickle(old_system.at, old_system.arrangePickle())
    print('Process aborted')
    print('Previous state succesfully saved')
    sys.exit()
except:
    out.savePickle(old_system.at, old_system.arrangePickle())
    print('Unexpected error')
    print('Previous state succesfully saved')
    raise
    pdb.set_trace()

out.savePickle(system.at, system.arrangePickle())
print('Simulation finished')
print('Last state succesfully saved')
sys.exit()

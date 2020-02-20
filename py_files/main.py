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
        self.at['part_solver'] = Leap_Frog(self.at['pic'])
        self.at['electrons'] = Electron_SW(0.0, c.E_SPWT, c.E_SIZE, c.DIM, c.DIM, self.at['mesh'].nPoints, c.NUM_TRACKED)
        self.at['protons'] = Proton_SW(0.0, c.P_SPWT, c.P_SIZE, c.DIM, c.DIM, self.at['mesh'].nPoints, c.NUM_TRACKED)
        self.at['e_field'] = Constant_Electric_Field(self.at['pic'], c.DIM)

    def arrangePickle(self):
        return ('ts', 'e_field', 'electrons', 'protons', 'part_solver')

    def arrangeVTK(self):
        return ('ts', 'e_field', 'electrons', 'protons')

#Initialization of the system and the previous step
system = System()
old_system = copy.deepcopy(system)

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
    pass

elif sys.argv[1] == '2':
    #File to be used as source of initial condition
    filename = 'ts1990.vtr'
    out.loadVTK(filename, system.at['mesh'], system.at, system.arrangeVTK())
    system.at['e_field'] = Electrostatic_2D_rm(system.at['pic'], c.DIM)

elif sys.argv[1] == '3':
    #File to be used as source of initial condition
    filename = 'sys_ts=57_2020-02-14_17h34m.pkl'
    out.loadPickle(filename, system.at, system.arrangePickle())

else:
    raise("Somehing is wrong here")


## ---------------------------------------------------------------------------------------------------------------
# Set up of the system before the Main loop
## ---------------------------------------------------------------------------------------------------------------

#Injection of particles at all outer boundaries
#Electrons
e_n = c.E_N*numpy.ones((len(system.at['pic'].mesh.boundaries[0].location)))
drift_e_vel = numpy.zeros((len(system.at['pic'].mesh.boundaries[0].location),2))
thermal_e_vel = c.E_V_TH*numpy.ones((len(system.at['pic'].mesh.boundaries[0].location)))
drift_e_vel[:,0] += c.E_V_SW
#Protons
p_n = c.P_N*numpy.ones((len(system.at['pic'].mesh.boundaries[0].location)))
drift_p_vel = numpy.zeros((len(system.at['pic'].mesh.boundaries[0].location),2))
thermal_p_vel = c.P_V_TH*numpy.ones((len(system.at['pic'].mesh.boundaries[0].location)))
drift_p_vel[:,0] += c.P_V_SW

for boundary in system.at['mesh'].boundaries:
    boundary.injectParticlesDummyBox(boundary.location, system.at['part_solver'], system.at['e_field'], system.at['electrons'], e_n, thermal_e_vel, drift_e_vel)
    boundary.injectParticlesDummyBox(boundary.location, system.at['part_solver'], system.at['e_field'], system.at['protons'], p_n, thermal_p_vel, drift_p_vel)

#Update of mesh values
system.at['part_solver'].updateMeshValues(system.at['protons'])
system.at['part_solver'].updateMeshValues(system.at['electrons'])

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
            system.at['e_field'].computeField([system.at['protons'], system.at['electrons']])
            system.at['part_solver'].advance(system.at['electrons'], system.at['e_field'])
            # Applying field borders and injecting at['electrons']
            for boundary in system.at['mesh'].boundaries:
                boundary.injectParticlesDummyBox(boundary.location, system.at['part_solver'], system.at['e_field'], system.at['electrons'], e_n, thermal_e_vel, drift_e_vel)
    
        #Proton motion
        system.at['part_solver'].advance(system.at['protons'], system.at['e_field'])
        for boundary in system.at['mesh'].boundaries:
            boundary.injectParticlesDummyBox(boundary.location, system.at['part_solver'], system.at['e_field'], system.at['protons'], p_n, thermal_p_vel, drift_p_vel)
    
        #Output vtk
        if system.at['ts']%10 == 0:
            out.saveVTK(system.at['mesh'], system.at, system.arrangeVTK())
        #if system.at['ts']%1 == 0:
        #    out.particleTracker(system.at['ts'], system.at['protons'], system.at['electrons'])
    
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

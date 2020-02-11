# Main file of the program
import copy
import numpy
import pdb
import sys
numpy.seterr(invalid='ignore', divide='ignore', over = 'raise')

import constants as c
from field import Constant_Electric_Field
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
#	+arrangePickle() : Variable = Return the attributes necessary to store the state of the system, in the order required.
#	+arrangeVTK() : Variable = Return the atrributes for save and load of/from VTK files, in the "order required".
class System(object):
    def __init__(self):
        self.ts = 0
        self.mesh = Mesh_2D_rm(c.XMIN, c.XMAX, c.YMIN, c.YMAX, c.NX, c.NY, c.DEPTH)
        self.pic = PIC_2D_rm1o(self.mesh)
        self.part_solver = Leap_Frog(self.pic)
        self.electrons = Electron_SW(0.0, c.E_SPWT, c.E_SIZE, c.DIM, c.DIM, self.mesh.nPoints, c.NUM_TRACKED)
        self.protons = Proton_SW(0.0, c.P_SPWT, c.P_SIZE, c.DIM, c.DIM, self.mesh.nPoints, c.NUM_TRACKED)
        self.e_field = Constant_Electric_Field(self.pic, c.DIM)

    def arrangePickle(self):
        return self.ts, self.e_field, self.electrons, self.protons, self.part_solver

    def arrangeVTK(self):
        return self.ts, self.e_field, self.electrons, self.protons

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
    filename = '.vtk'
    out.loadVTK(filename, system.mesh, *system.arrangeVTK())

elif sys.argv[1] == '3':
    #File to be used as source of initial condition
    filename = '.pkl'
    out.loadPickle(filename, *arrangePickle())

else:
    raise("Somehing is wrong here")


## ---------------------------------------------------------------------------------------------------------------
# Set up of the system before the Main loop
## ---------------------------------------------------------------------------------------------------------------

#Injection of particles at all outer boundaries
#Electrons
e_n = c.E_N*numpy.ones((len(system.pic.mesh.boundaries[0].location)))
drift_e_vel = numpy.zeros((len(system.pic.mesh.boundaries[0].location),2))
thermal_e_vel = c.E_V_TH*numpy.ones((len(system.pic.mesh.boundaries[0].location)))
drift_e_vel[:,0] += c.E_V_SW
#Protons
p_n = c.P_N*numpy.ones((len(system.pic.mesh.boundaries[0].location)))
drift_p_vel = numpy.zeros((len(system.pic.mesh.boundaries[0].location),2))
thermal_p_vel = c.P_V_TH*numpy.ones((len(system.pic.mesh.boundaries[0].location)))
drift_p_vel[:,0] += c.P_V_SW

for boundary in system.mesh.boundaries:
    boundary.injectParticlesDummyBox(boundary.location, system.part_solver, system.e_field, system.electrons, e_n, thermal_e_vel, drift_e_vel)
    boundary.injectParticlesDummyBox(boundary.location, system.part_solver, system.e_field, system.protons, p_n, thermal_p_vel, drift_p_vel)

#Half-step speed recoil to meet Leap-frog condition
system.part_solver.updateMeshValues(system.protons)
system.part_solver.updateMeshValues(system.electrons)

## ---------------------------------------------------------------------------------------------------------------
# Main loop
## ---------------------------------------------------------------------------------------------------------------

# try/except block to capture any error and save before the step before the crash. It also allows the simulation to be stopped.
try:
    while system.ts < c.NUM_TS:
        print('ts = ', system.ts)
    
        # Electron motion
        for te in range(c.ELECTRON_TS):
            print('te = ', te)
            system.e_field.computeField([system.protons, system.electrons])
            system.part_solver.advance(system.electrons, system.e_field)
            # Applying field borders and injecting electrons
            for boundary in system.mesh.boundaries:
                boundary.injectParticlesDummyBox(boundary.location, system.part_solver, system.e_field, system.electrons, e_n, thermal_e_vel, drift_e_vel)
    
        #Proton motion
        system.part_solver.advance(system.protons, system.e_field)
        for boundary in system.mesh.boundaries:
            boundary.injectParticlesDummyBox(boundary.location, system.part_solver, system.e_field, system.protons, p_n, thermal_p_vel, drift_p_vel)
    
        #Output vtk
        if system.ts%10 == 0:
            out.saveVTK(system.mesh, *system.arrangeVTK())
        if system.ts%1 == 0:
            out.particleTracker(system.ts, system.protons, system.electrons)
    
        #Updating previous state
        old_system = copy.deepcopy(system)
    
        #Advance in timestep
        system.ts += 1

except KeyboardInterrupt:
    out.savePickle(*old_system.arrangePickle())
    print('Process aborted')
    print('Previous state succesfully saved')
    sys.exit()
except:
    out.savePickle(*old_system.arrangePickle())
    print('Unexpected error')
    print('Previous state succesfully saved')
    raise
    pdb.set_trace()

out.savePickle(*system.arrangePickle())
print('Simulation finished')
print('Last state succesfully saved')
sys.exit()

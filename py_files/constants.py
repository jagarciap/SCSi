import numpy

#Physical constants
EPS_0 =8.85418782e-12   # F/m, vacuum permittivity
K     =1.38064825e-23          # J/K, Boltzmann constant
ME    =9.10938215e-31       # kg, electron mass
QE    = -1.602176565e-19      # C, electron charge
AMU   =1.660538921e-27  # kg, atomic mass unit
MP = 1.007276467*AMU # kg, proton mass
EV_TO_K=11604.525        # 1eV in Kelvin
g = 9.80665 #Gravity on Earth surface

#Simulation time's parameters
NUM_TS = numpy.int(4e3)   #Number of total steps in the system
P_DT= 5e-8            # time step size
E_DT = 5e-9        # time step for electron dynamics
ELECTRON_TS = 10   #Internal iterations for electron dynamics

#Geometrical parameters for a rectangular outer boundary
XMIN = 0.0
XMAX = 2.0#5.0
YMIN = -1.0#-205
YMAX = 1.0#2.5
DEPTH = 1.0
NX = numpy.uint16(41)#10
NY = numpy.uint16(41)#10

#Particle physical parameters
E_N = 7e9
E_T = 84.47*EV_TO_K
E_V_TH = numpy.sqrt(2*K*E_T/ME)
E_V_SW = 300e3

P_N = 7e9
P_T = 87.25*EV_TO_K
P_V_TH = numpy.sqrt(2*K*P_T/MP)
P_V_SW = 300e3

#Particle simulation parameters
P_SIZE = numpy.uint32(1e5)     #Size of the ions array
E_SIZE = numpy.uint32(1e5)     #Size of the electrons array
P_SPWT = 4e5#5e6
E_SPWT = 4e5#5e6

#Number of particles traced
NUM_TRACKED = 100
DIM = numpy.uint8(2)

import numpy as np

#Physical constants
EPS_0 =8.85418782e-12   # F/m, vacuum permittivity
K     =1.38064825e-23          # J/K, Boltzmann constant
ME    =9.10938215e-31       # kg, electron mass
QE    = -1.602176565e-19      # C, electron charge
AMU   =1.660538921e-27  # kg, atomic mass unit
EV_TO_K=11604.525        # 1eV in Kelvin
g = 9.80665 #Gravity on Earth surface

#Simulation time's parameters
NUM_TS = 1000   #Number of total steps in the system
ELECTRON_TS = 250   #Internal iterations for electron dynamics
DT= 5e-8            # time step size
E_DT = 2e-10        # time step for electron dynamics

#Geometrical parameters
R_MIN = -5e-2#3.5e-2        # r min para SPT
R_MAX = 5e-2        # r max para SPT
Z_MIN = -5e-2             # z min para SPT
Z_MAX = 5e-2#2.5E-2      # acceleration channel exit
NR = 40#12             # cell number in r-direction
NZ = 40#50             # cell number in z-direction
DR= (R_MAX-R_MIN)/NR         # cell r-spacing
DZ= (Z_MAX-Z_MIN)/NZ#6.25e-4         # cell z-spacing

#Particle parameters
I_SIZE = 1     #Size of the ions array
E_SIZE = 1     #Size of the electrons array
I_SPWT = 1
E_SPWT = 1

#Number of particles traced
NUM_PART = E_SIZE

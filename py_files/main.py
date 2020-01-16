# Main file of the program
import numpy
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

mesh = Mesh_2D_rm(c.XMIN, c.XMAX, c.YMIN, c.YMAX, c.NX, c.NY, c.DEPTH)
pic = PIC_2D_rm1o(mesh)
e_field = Constant_Electric_Field(pic, mesh.boundaries, mesh.nPoints, c.DIM)
electrons = Electron_SW(0.0, c.E_SPWT, c.E_SIZE, c.DIM, c.DIM, mesh.nPoints)
protons = Proton_SW(0.0, c.P_SPWT, c.P_SIZE, c.DIM, c.DIM, mesh.nPoints)
part_solver = Leap_Frog(pic)

## ---------------------------------------------------------------------------------------------------------------
# Initial condition
## ---------------------------------------------------------------------------------------------------------------

#Injection of particles in the left side and conditions for velocity
left_boundary = numpy.arange(0,c.NX*c.NY,c.NX)
#Electrons
drift_e_vel = numpy.zeros((len(left_boundary),2))
thermal_e_vel = c.E_V_TH*numpy.ones((len(left_boundary)))
drift_e_vel[:,0] += c.E_V_SW
#Protons
drift_p_vel = numpy.zeros((len(left_boundary),2))
thermal_p_vel = c.P_V_TH*numpy.ones((len(left_boundary)))
drift_p_vel[:,0] += c.P_V_SW

for boundary in mesh.boundaries:
    boundary.applyElectricBoundary(e_field)
    boundary.injectParticlesDummyBox(pic.mesh.boundaries[0].location, pic, electrons, c.E_N, thermal_e_vel, drift_e_vel)
    boundary.injectParticlesDummyBox(pic.mesh.boundaries[0].location, pic, protons, c.P_N, thermal_p_vel, drift_p_vel)
part_solver.initialConfiguration(protons, e_field)
part_solver.initialConfiguration(electrons, e_field)

## ---------------------------------------------------------------------------------------------------------------
# Main loop
## ---------------------------------------------------------------------------------------------------------------

for tp in range(c.NUM_TS):
    print('tp = ', tp)

    # Electron motion
    for te in range(c.ELECTRON_TS):
        print('te = ', te)
        e_field.computeField([protons, electrons])
        part_solver.advance(electrons, e_field)
        # Applying field borders and injecting electrons
        for boundary in mesh.boundaries:
            boundary.applyElectricBoundary(e_field)
            boundary.injectParticlesDummyBox(pic.mesh.boundaries[0].location, pic, electrons, c.E_N, thermal_e_vel, drift_e_vel)

    #Proton motion
    part_solver.advance(protons, e_field)
    for boundary in mesh.boundaries:
        boundary.injectParticlesDummyBox(pic.mesh.boundaries[0].location, pic, protons, c.P_N, thermal_p_vel, drift_p_vel)

    #Output
    if tp%10 == 0:
        out.output_vtk(tp, mesh, electrons, protons, e_field)

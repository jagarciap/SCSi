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
e_field = Constant_Electric_Field(pic, c.DIM)
electrons = Electron_SW(0.0, c.E_SPWT, c.E_SIZE, c.DIM, c.DIM, mesh.nPoints, c.NUM_TRACKED)
protons = Proton_SW(0.0, c.P_SPWT, c.P_SIZE, c.DIM, c.DIM, mesh.nPoints, c.NUM_TRACKED)
part_solver = Leap_Frog(pic)

## ---------------------------------------------------------------------------------------------------------------
# Initial condition
## ---------------------------------------------------------------------------------------------------------------

#Injection of particles at all outer boundaries
#Electrons
e_n = c.E_N*numpy.ones((len(pic.mesh.boundaries[0].location)))
drift_e_vel = numpy.zeros((len(pic.mesh.boundaries[0].location),2))
thermal_e_vel = c.E_V_TH*numpy.ones((len(pic.mesh.boundaries[0].location)))
drift_e_vel[:,0] += c.E_V_SW
#Protons
p_n = c.P_N*numpy.ones((len(pic.mesh.boundaries[0].location)))
drift_p_vel = numpy.zeros((len(pic.mesh.boundaries[0].location),2))
thermal_p_vel = c.P_V_TH*numpy.ones((len(pic.mesh.boundaries[0].location)))
drift_p_vel[:,0] += c.P_V_SW

for boundary in mesh.boundaries:
    boundary.injectParticlesDummyBox(boundary.location, pic, electrons, e_n, thermal_e_vel, drift_e_vel)
    boundary.injectParticlesDummyBox(boundary.location, pic, protons, p_n, thermal_p_vel, drift_p_vel)
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
            boundary.injectParticlesDummyBox(boundary.location, pic, electrons, e_n, thermal_e_vel, drift_e_vel)

    #Proton motion
    part_solver.advance(protons, e_field)
    for boundary in mesh.boundaries:
        boundary.injectParticlesDummyBox(boundary.location, pic, protons, p_n, thermal_p_vel, drift_p_vel)

    #Output vtk
    if tp%10 == 0:
        out.output_vtk(tp, mesh, electrons, protons, e_field)
    if tp%1 == 0:
        out.particle_tracker(tp, protons, electrons)

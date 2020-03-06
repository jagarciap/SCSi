# Data structures for electrons
from Species.species import Species
import constants as c


#Electron (Inherits from Species):
#
#Definition = Species that take care of electrons.
#Attributes:
#	+Species attributes.
class Electron(Species):
    def __init__(self, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints, n_num_tracked, n_string):
        super().__init__("Electron"+n_string, c.E_DT ,c.QE, c.ME, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints, n_num_tracked)

#Electron_SW (Inherits from Species):
#
#Definition = Species that take care of electrons coming from solar wind.
#Attributes:
#	+type (string) = "Electron - Solar wind"
#	+Species attributes.
class Electron_SW(Electron):
    def __init__(self, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints, n_num_tracked = 0):
        super().__init__(n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints, n_num_tracked, " - Solar wind")

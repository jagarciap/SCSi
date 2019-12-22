# Data structures for electrons
from Species.species import Species
import constants as c


#Electron (Inherits from Species):
#
#Definition = Species that take care of electrons.
#Attributes:
#	+type (string) = some string descriptor that indicate the type/source of electrons.
#	+Species attributes.
class Electron(Species):
    def __init__(self, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints, n_string):
        self.type = n_string
        super().__init__(c.E_DT ,c.QE, c.ME, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints)

#Electron_SW (Inherits from Species):
#
#Definition = Species that take care of electrons coming from solar wind.
#Attributes:
#	+type (string) = "Electron - Solar wind"
#	+Species attributes.
class Electron_SW(Electron):
    def __init__(self, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints):
        super().__init__(n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints, "Electron - Solar wind")

# Data structures for protons
from Species.species import Species
import constants as c

#Proton (Inherits from Species):
#
#Definition = Species that take care of protons.
#Attributes:
#	+type (string) = some string descriptor that indicate the type/source of protons.
#	+Species attributes.
class Proton(Species):
    def __init__(self, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints, n_string):
        self.type = n_string
        super().__init__(c.P_DT, -c.QE, c.MP, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints)

#Proton_SW (Inherits from Species):
#
#Definition = Species that take care of protons coming from solar wind.
#Attributes:
#	+type (string) = "Proton - Solar wind"
#	+Species attributes.
class Proton_SW(Proton):
    def __init__(self, n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints):
        super().__init__(n_debye, n_spwt, n_max_n, n_pos_dim, n_vel_dim, n_nPoints, "Proton - Solar wind")

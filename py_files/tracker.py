#File that handles the classes for tracking particles
import numpy

#Tracker:
#
#Defintion = Class that initializes and stores the information for the tool of tracking particles after execution.
#Attributes:
#	+trackers ([Species_Tracker]) = List of instances of Species_Tracker.
#	+num_tracked (int) = Number of particles to be tracked.
#Methods:
#	+print() = Creates the files and print the information of the particles. Is basically the function to be called for output of information.
class Tracker(object):
    def __init__(self, num_tracked, *args):
        self.trackers = []
        self.num_tracked = numpy.uint8(num_tracked)
        for species in args:
            species_tracker = Species_Tracker(species, num_tracked)
            self.trackers.append(species_tracker)

    def print(self):
        # Checking tracking method
        for spc in self.trackers:
            if numpy.any(numpy.isnan(spc.indices)):
                raise ValueError("There should not be any nan values")

        ind_ions = numpy.flatnonzero(ions.part.printable[:ions.np])
        ind_neutrals = numpy.flatnonzero(neutrals.part.printable[:neutrals.np])
        np_ions = len(ind_ions)
        np_neutrals = len(ind_neutrals)
        if len(ions.part.print_ind) > 0 or len(neutrals.part.print_ind) > 0:
            pdb.set_trace()

        narray = numpy.zeros((self.num_tracked,2*len(self.trackers)))
        narray[ions.part.printable[ind_ions]-1, 0] = ions.part.x[ind_ions,0]
        narray[ions.part.printable[ind_ions]-1, 1] = ions.part.x[ind_ions,1]
        narray[neutrals.part.printable[ind_neutrals]-1, 2] = neutrals.part.x[ind_neutrals,0]
        narray[neutrals.part.printable[ind_neutrals]-1, 3] = neutrals.part.x[ind_neutrals,1]

        cwd = os.getcwd()
        workfile = cwd+'/particle_tracker/ts={:05d}.dat'.format(ts)
        nheader = 'Max. number of particles: \t {:4d} \n Ions \t Neutrals \n {:3d} \t {:3d} \n'.format( num,np_ions, np_neutrals)
        numpy.savetxt(workfile, narray , fmt = '%.5e', delimiter = '\t', header = nheader)

#Species_Tracker (Composition with Tracker):
#
#Definition = Class that store the indices of the particles being tracked for a certain species.
#Attributes:
#	+identifier (String) = Same species.type
#	+indices (int) = array of size num_tracked that store the indices of the particles as stored in Particles class.
class Species_Tracker(object):
    def __init__(self, species, num_tracked):
        self.identifier = species.type
        indices = numpy.nan*numpy.ones((num_tracked), dtype = numpy.uint32)

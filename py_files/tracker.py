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


        narray = numpy.zeros((self.num_tracked, numpy.shape(self.trackers[0].position)[1]*len(self.trackers)))
        nHeader = ''
        for i in len(self.trackers):
            narray[:,2*i:2*(i+1)] = self.trackers[i].position[self.trackers[i].indices,:]
            nHeader += self.trackers[i].identifier + '\t'

        cwd = os.path.split(os.getcwd())[0]
        workfile = cwd+'/particle_tracker/ts={:05d}.dat'.format(ts)
        nHeader = 'No. of particles = {:d} \n'.format(self.num_tracked)+nHeader+'\n' 
        numpy.savetxt(workfile, narray , fmt = '%.5e', delimiter = '\t', header = nHeader)

#Species_Tracker (Composition with Tracker):
#
#Definition = Class that store the indices of the particles being tracked for a certain species.
#Attributes:
#	+identifier (String) = Same species.type
#	+indices (int) = array of size num_tracked that store the indices of the particles as stored in Particles class.
#       +position ([double,double]) = reference to the list of positions of the species being tracked.
class Species_Tracker(object):
    def __init__(self, species, num_tracked):
        self.identifier = species.type
        self.indices = numpy.nan*numpy.ones((num_tracked), dtype = numpy.uint32)
        self.position = species.part_values.position

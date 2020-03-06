

#---------------------------------------------------------------------------------------------------------------------------------------
# #1: I was trying to write the function using numpy instead of for
#---------------------------------------------------------------------------------------------------------------------------------------
#       +createDummyBox([ind]location, PIC pic, Species species, [double] delta_n, [double] n_vel, [double] shift_vel) = create the dummy boxes with particles in them.
    def createDummyBox(self, location, pic, species, delta_n, n_vel, shift_vel):
        # Locating borders
        phys_loc = pic.mesh.getPosition(location)
        left = numpy.equal(phys_loc[:,0], pic.mesh.xmin)
        right = numpy.equal(phys_loc[:,0], pic.mesh.xmax)
        bottom = numpy.equal(phys_loc[:,1], pic.mesh.ymin)
        top = numpy.equal(phys_loc[:,1], pic.mesh.ymax)
        # Corners
        lb = numpy.logical_and(left, bottom) 
        lt = numpy.logical_and(left, top)
        rb = numpy.logical_and(right, bottom)
        rt = numpy.logical_and(right, top)
        # Margins
        corners = numpy.logical_or(lb, numpy.logical_or(lt, numpy.logical_or(rb, rt)))
        nc = numpy.logical_not(corners)
        l = numpy.logical_and(left, nc)
        r = numpy.logical_and(right, nc)
        t = numpy.logical_and(top, nc)
        b = numpy.logical_and(bottom,nc)
        # Creation of particles (using the short version of adding randomness)
        mpf_new = delta_n*pic.mesh.volumes[location]/species.spwt
        mpf_new = numpy.where(corners, mp_new*3)
        mp_new = mpf_new+numpy.random.rand(len(location))).astype(int)
        # Generate this many particles
        total_new = numpy.sum(mp_new)
        pos = numpy.zeros((total_new, species.pos_dim))
        vel = numpy.zeros((total_new, species.vel_dim))
        # Setting up velocities and positions
        c = 0
        for i in range(len(location)):
            vel[c:c+mp_new[i],:] += self.sampleIsotropicVelocity(n_vel[i], mp_new[i])+shift_vel[i,:]
            if phys_loc[i,0] == mesh.xmin:
                if phys_loc[i,1] == mesh.ymin:
                    pos[c:c+mp_new[i],0] += phys_loc[i,0] + numpy.random.rand(mp_new[i])*box_x[i]
            #pos[c:c+mp_new[i],0] += phys_loc[i,0] + (numpy.random.rand(mp_new[i])-0.5)*pic.mesh.dx
            pos[c:c+mp_new[i],1] += phys_loc[i,1] + (numpy.random.rand(mp_new[i])-0.5)*pic.mesh.dy
            c += mp_new[i]

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

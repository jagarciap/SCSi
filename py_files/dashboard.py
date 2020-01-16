

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

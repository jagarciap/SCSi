#File containing the methods to calculate potentials and fields
import numpy
import pdb
#NOTE: Delete this later
import matplotlib.pyplot as plt


#       +Method that computes the electric potential for a 2D rectangular mesh using the method of
#       +Successive Over Relaxation with Chebyshev Acceleration (SORCA). The method assumes that the 
#       +boundaries has been treated previously and assumes a Dirichlet boundary, so it only calculates values
#       +for the inner nodes. If a Neumann condition is desired, then the outer layer nodes are calculated correspondingly 
#       +in the later method apply applyElectricBoundary(e_field), not here.
#       +Parameters:
#       +rhs ([double]) =  is the right-hand-side of the equation Poisson(\phi) = rhs.
#       +pot ([double]) = is the potential before being updated, or it can be seen as the first guess of the solution.
#       +err (duble) = is the maximum variation expected between succesives steps to reach solution. 
#       +step_limit (int) = a boundary in number of steps to avoid infinite loops.
def poissonSolver_2D_rm_SORCA(mesh, pot, rhs, err = 1e-3, step_limit = 1000):
    #Storing 1D to 2D conversion
    ind_2D = mesh.arrayToIndex(numpy.arange(mesh.nPoints))
    #Coefficients of the equation
    a = 1/mesh.dy/mesh.dy*numpy.ones(mesh.nPoints)
    b = a
    c = 1/mesh.dx/mesh.dx*numpy.ones(mesh.nPoints)
    d = c
    e = -2*(1/mesh.dx/mesh.dx+1/mesh.dy/mesh.dy)*numpy.ones(mesh.nPoints)
    #Defining rho
    delta_sq = mesh.nx*mesh.nx if mesh.nx < mesh.ny else mesh.ny*mesh.ny
    rho = 1 - numpy.pi*numpy.pi/2/delta_sq
    #NOTE: Delete later
    tot_norm = []
    #Solver
    w = 1.0
    #w = 2/(1-numpy.pi/mesh.nx)
    for t in range(1, step_limit+1):
        norm = 0.0
        for ind in range (mesh.nPoints):
            bound = False
            for boundary in mesh.boundaries:
                if ind in boundary.location:
                    bound = True
                    break
            if (ind%mesh.nx+ind//mesh.nx)%2 == t%2 or bound:
                continue
            res = a[ind]*pot[ind-mesh.nx]+b[ind]*pot[ind+mesh.nx]+c[ind]*pot[ind-1]+d[ind]*pot[ind+1]+e[ind]*pot[ind]-rhs[ind]
            norm += res*res
            pot[ind] = pot[ind]-w*res/e[ind]
        tot_norm.append(norm)
        #print("norm",norm, "w", w)
        #mid = int(mesh.ny/2)
        #ext = int(mesh.nx/4)
        #for j in range (-ext+mid, mid+ext+1):
        #    print(pot[j*mesh.nx+mid-ext:j*mesh.nx+mid+ext+1])
        #print("-------------------------------")
        #pdb.set_trace()
        w = 1.0/(1-0.25*rho*rho*w)
        if t == 1:
            w = 1.0/(1-0.5*rho*rho)
        if norm < err:
            break
        if t == step_limit:
            pdb.set_trace()
            raise ValueError("step_limit reached, solution not obtained. Error = {:e}.".format(norm))

#       +Derivation of the scalar field potential ([double]) with the method of central differences, at nodes not in the boundaries.
def derive_2D_rm(mesh, potential):
    #Creating array of indexes to slice through
    ind = numpy.arange(mesh.nPoints, dtype = numpy.uint16)
    for boundary in mesh.boundaries:
        ind = numpy.delete( ind, boundary.location)
    #Creating temporary field
    field = numpy.zeros((mesh.nPoints, 2))
    #Derivation
    field[ind,0] = (potential[ind+1]-potential[ind-1])/(2*mesh.dx)
    field[ind,1] = (potential[ind+mesh.nx]-potential[ind-mesh.nx])/(2*mesh.dy)
    return field

#       +Pade derivation for the nodes in the boundaries. Normal 2nd order derivation when the boundary is perpendicular to the direction of the derivative.
#       +Arguments:
#       +location ([ind]) = location of the boundary nodes to be treated.
#       +mesh (Outer_2D_Rectangular) = mesh with the information to make the finite difference.
#       +potential ([double]) = scalar to be derivated.
#       +Return: [double, double] two-component derivation of potential, with every row being one node of location.
def derive_2D_rm_boundaries(location, mesh, potential):
    #Creating temporary field
    field = numpy.zeros((len(location),2))
    for ind in range(len(location)):
        #Handling corners
        if location[ind] == 0:
            field[ind,0] = (-3*potential[location[ind]]+4*potential[location[ind]+1]-potential[location[ind]+2])/(2*mesh.dx)
            field[ind,1] = (-3*potential[location[ind]]+4*potential[location[ind]+mesh.nx]-potential[location[ind]+2*mesh.nx])/(2*mesh.dy)
        elif location[ind] == mesh.nx:
            field[ind,0] = (3*potential[location[ind]]-4*potential[location[ind]-1]+potential[location[ind]-2])/(2*mesh.dx)
            field[ind,1] = (-3*potential[location[ind]]+4*potential[location[ind]+mesh.nx]-potential[location[ind]+2*mesh.nx])/(2*mesh.dy)
        elif location[ind] == mesh.nx*(mesh.ny-1):
            field[ind,0] = (-3*potential[location[ind]]+4*potential[location[ind]+1]-potential[location[ind]+2])/(2*mesh.dx)
            field[ind,1] = (3*potential[location[ind]]-4*potential[location[ind]-mesh.nx]+potential[location[ind]-2*mesh.nx])/(2*mesh.dy)
        elif location[ind] == mesh.nx*mesh.ny-1:
            field[ind,0] = (3*potential[location[ind]]-4*potential[location[ind]-1]+potential[location[ind]-2])/(2*mesh.dx)
            field[ind,1] = (3*potential[location[ind]]-4*potential[location[ind]-mesh.nx]+potential[location[ind]-2*mesh.nx])/(2*mesh.dy)
        #Handling non-corner borders
        elif location[ind] < mesh.nx:
            field[ind,0] = (potential[location[ind]+1]-potential[location[ind]-1])/(2*mesh.dx)
            field[ind,1] = (-3*potential[location[ind]]+4*potential[location[ind]+mesh.nx]-potential[location[ind]+2*mesh.nx])/(2*mesh.dy)
        elif location[ind]%mesh.nx == 0:
            field[ind,0] = (-3*potential[location[ind]]+4*potential[location[ind]+1]-potential[location[ind]+2])/(2*mesh.dx)
            field[ind,1] = (potential[location[ind]+mesh.nx]-potential[location[ind]-mesh.nx])/(2*mesh.dy)
        elif location[ind]%mesh.nx == mesh.nx-1:
            field[ind,0] = (3*potential[location[ind]]-4*potential[location[ind]-1]+potential[location[ind]-2])/(2*mesh.dx)
            field[ind,1] = (potential[location[ind]+mesh.nx]-potential[location[ind]-mesh.nx])/(2*mesh.dy)
        else:
            field[ind,0] = (potential[location[ind]+1]-potential[location[ind]-1])/(2*mesh.dx)
            field[ind,1] = (3*potential[location[ind]]-4*potential[location[ind]-mesh.nx]+potential[location[ind]-2*mesh.nx])/(2*mesh.dy)
    return field

#File containing the methods to calculate potentials and fields
import numpy


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
    delta_sq = dx*dx if dx < dy else dy*dy
    rho = 1 - numpy.pi*numpy.pi/delta_sq
    #Solver
    w = 1.0
    for t in range(1, step_limit+1):
        norm = 0.0
        for ind in range (mesh.nPoints):
            bound = False
            for boundary in mesh.boundaries:
                if ind in boundary.location:
                    bound = True
                    break
            if ind%2 == t%2 or bound:
                continue
            res = a[ind]*pot[ind-mesh.nx]+b[ind]*pot[ind+mesh.nx]+c[ind]*pot[ind-1]+d[ind]*pot[ind+1]+e[ind]*pot[ind]-rhs[ind]
            norm += res*res
            pot[ind] = pot[ind]-w*res/e[ind]
        w = 1.0/(1-0.25*rho*rho*w)
        if t == 1:
            w = 1.0/(1-0.5*rho*rho)
        if norm < err:
            break
        if t == step_limit:
            raise ValueError("step_limit reached, solution not obtained. Error = {:e}.".format(norm))


def derive_2D_rm(mesh, self.potential):


import matplotlib.pyplot as plt
import numpy as np
import math
import time
import meshio
import sys
from scipy.spatial import distance_matrix as dist_mat
from scipy import linalg
from scipy.spatial import Delaunay
from dolfin import *

max_dist = np.linalg.norm([141.7/2, 141.7/2, 60.5])
halves = [141.7/2, 141.7/2, 60]
lower_bound = .01
upper_bound = 100
b_slope = (upper_bound - lower_bound)/max_dist
A = (upper_bound - lower_bound)/(max_dist**2)
sphere_r = 10

class inside_surface(SubDomain):
    def inside(self, x, on_boundary):
        return np.linalg.norm(x - halves) <= sphere_r + .1

class mu_nw(UserExpression):
    def __init__(self, mesh, nu, cell_points, **kwargs):
        self.mesh = mesh
        self.nu = nu
        self.cell_points = cell_points
        self._mins = np.ones(3)*100
        self._maxs = np.zeros(3)
        super().__init__(**kwargs)

    def eval(self, value, x):
        dist2cell = np.linalg.norm(self.cell_points - x, axis=1).min()
        # _E = A*(dist2cell**2) + lower_bound
        _E = b_slope * dist2cell + lower_bound
        # _E = 100


        value[0] = _E / (2 * (1 + self.nu))

class lmbda_nw(UserExpression):
    def __init__(self, mesh, nu, cell_points, **kwargs):
        self.mesh = mesh
        self.nu = nu
        self.cell_points = cell_points
        super().__init__(**kwargs)

    def eval(self, value, x):
        dist2cell = np.linalg.norm(self.cell_points - x, axis=1).min()
        # _E = A*(dist2cell**2) + lower_bound
        _E = b_slope * dist2cell + lower_bound
        # _E = 100

        value[0] = _E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))

def write_interp2mesh(mesh, interpolation_PETSc, file_in, save_tag):
    print("writing interp2mesh")
    coordinates = mesh.coordinates()
    cells = mesh.cells()
    interpolation = np.asarray([interpolation_PETSc.vector()[ii] for ii in range(len(interpolation_PETSc.vector()))])
    centers = np.asarray([coordinates[cell].mean(0) for cell in cells])
    val_per_cell = np.asarray([interpolation_PETSc(cent) for cent in centers])
    mesh = meshio.Mesh(points=coordinates,
                       cells={"tetra": cells},
                       cell_data={"tetra": {"geometrical": val_per_cell}})
    meshio.write("saved_interp/" + file_in + "_" + save_tag + ".vtk", mesh)
    np.savetxt("saved_interp/" + file_in + "_" + save_tag + "_points.txt", coordinates, delimiter=",")
    np.savetxt("saved_interp/" + file_in + "_" + save_tag + "_mu.txt", interpolation, delimiter=",")
    print("written interp2mesh")
    return interpolation, centers, val_per_cell

file_in = sys.argv[1]
msh_file = sys.argv[2]
shrink_r = float(sys.argv[3])
save_tag = sys.argv[4]

mod_const_flag = True
xy_res = 141.7/512
z_res = 120/151
xy_bound = 141.7
z_bound = 120
nbeads = 4000

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"

cell_points = meshio.read("MeshesMSH/" + msh_file + ".msh").points
mesh = Mesh()
with XDMFFile("MeshesXDMF/" + file_in + ".xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 3)
with XDMFFile("MeshesXDMF/" + file_in + "_function.xdmf") as infile:
    infile.read(mvc, "geometrical")
subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvc)

gmark, cmark = 0, 100
dx = Measure('dx',domain=mesh, subdomain_data=subdomains, metadata={'quadrature_degree': 2})

V = VectorFunctionSpace(mesh, "Lagrange", 1)
sF = FunctionSpace(mesh, "Lagrange", 1)
du = TrialFunction(V)
v  = TestFunction(V)
u  = Function(V)

inside = inside_surface()
inside.mark(subdomains, 500)

zero = Constant((0.0, 0.0, 0.0))
bcs = []
sbd = []
sbd.append(CompiledSubDomain("near(x[0], side)", side = 0))
sbd.append(CompiledSubDomain("near(x[1], side)", side = 0))
sbd.append(CompiledSubDomain("near(x[2], side)", side = 0))
sbd.append(CompiledSubDomain("near(x[0], side)", side = xy_bound))
sbd.append(CompiledSubDomain("near(x[1], side)", side = xy_bound))
sbd.append(CompiledSubDomain("near(x[2], side)", side = z_bound))
[bcs.append((DirichletBC(V, zero, sub))) for sub in sbd]

shrink_list = [.95, .9, .85, .8, .75, .7]

shr0, nu0 = 1.0, 0.45

mu_inst = mu_nw(mesh, nu0, cell_points)
mu = interpolate(mu_inst, sF)

lmbda_inst = lmbda_nw(mesh, nu0, cell_points)
lmbda = interpolate(lmbda_inst, sF)

interpolation, centers, val_per_cell = write_interp2mesh(mesh, mu, file_in, save_tag)

r = Expression(("shrink_r*((halves_x - x[0])/sphere_r)",
                "shrink_r*((halves_y - x[1])/sphere_r)",
                "shrink_r*((halves_z - x[2])/sphere_r)"),
               shrink_r=shrink_r,
               halves_x=halves[0],
               halves_y=halves[1],
               halves_z=halves[2],
               sphere_r=sphere_r,
               degree=2)

bcs.append(DirichletBC(V, r, inside))

d = len(u)
I = Identity(d)
F = I + grad(u)
B = Constant((0.0, 0.0, 0.0))
T = Constant((0.0, 0.0, 0.0))
C = F.T * F

# if mod_const:
#     E0 = 1
#     mu0 = Constant(E0 / (2 * (1 + nu0)))
#     lmbda0 = Constant(E0 * nu0 / ((1 + nu0) * (1 - 2 * nu0)))
# else:
#     mu0 = mu_nw(mesh, nu0, cell_points)
#     lmbda0 = lmbda_nw(mesh, nu0, cell_points)


Ic = tr(C)
J = det(F)

psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J)) ** 2
Pi = psi * dx - dot(B, u) * dx - dot(T, u) * ds
F = derivative(Pi, u, v)
J = derivative(F, u, du)

solve_start = time.time()
solve(F == 0, u, bcs, J=J, solver_parameters={"newton_solver": {"linear_solver": "mumps"}})
print("Solve time: ", time.time() - solve_start)

hdf5_file = HDF5File(mesh.mpi_comm(),
                     "saved_functions/" + file_in + "_" + str(shrink_r) + "_" + save_tag + "_func.h5",
                     "w")
hdf5_file.write(u, "/function")
hdf5_file.close()

file = File("Solutions/" + file_in + "_" + str(shrink_r) + "_" + save_tag + ".pvd")
file << u


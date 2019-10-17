import matplotlib.pyplot as plt
import numpy as np
import time
import meshio
import sys
from scipy import linalg
from dolfin import *


def sim_beads(mesh, u, subdomains, gmark, n, xy_bound, z_bound):
    sub_arr = subdomains.array()
    mcells = mesh.cells()
    n_cells, _ = mcells.shape

    bounds = np.array([[.1, .1, .1],
                       [xy_bound, xy_bound, z_bound]])
    beads = []
    displacements = []

    valid_count = 0
    while valid_count<n:
        bead_candidate = np.random.uniform(bounds[0], bounds[1], size=(3))
        candidate_idx = mesh.bounding_box_tree().compute_first_entity_collision(Point(bead_candidate[0],
                                                                                    bead_candidate[1],
                                                                                    bead_candidate[2]))
        if sub_arr[candidate_idx] == gmark:
            displacements.append(u(bead_candidate))
            valid_count = valid_count + 1
            beads.append(bead_candidate)

    return np.array(beads), np.array(displacements)

file_in = sys.argv[1]
save_tag = sys.argv[2]
mod_const_flag = False
xy_res = 141.7/512
z_res = 120/151
xy_bound = 141.7
z_bound = 120
nbeads = 4000

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
du = TrialFunction(V)
v  = TestFunction(V)
u  = Function(V)

shrink_list = [.9, .85, .8, .75, .7]

start = time.time()
for shrink_param in shrink_list:

    hdf5_file = HDF5File(mesh.mpi_comm(),
                         "saved_functions/" + file_in + "_" + str(shrink_param).ljust(4, '0') + "_" + save_tag + "_func.h5",
                         "r")
    hdf5_file.read(u, "/function")
    hdf5_file.close()
    beads, displacements = sim_beads(mesh, u, subdomains, gmark, nbeads, xy_bound, z_bound)
    np.savetxt("beads_from_func/" + file_in + "_" + str(shrink_param).ljust(4, '0') + "_" + save_tag + "_beads.txt", beads, delimiter=",")
    np.savetxt("beads_from_func/" + file_in + "_" + str(shrink_param).ljust(4, '0') + "_" + save_tag + "_disps.txt", displacements, delimiter=",")



import numpy as np
import matplotlib.pyplot as plt
import sys
from dolfin import *


def tags2labels(tags):
    label_list = []
    for tag in tags:
        split_out = tag.split('_')
        if len(split_out) == 3:
            if 'old' in tag:
                if split_out[-1] == '1111':
                    label_list.append("E = 0.1 Old Model")
                else:
                    label_list.append("E = " + split_out[-1] + " Old Model")
            else:
                if split_out[-1] == '1111':
                    label_list.append("E = 0.1")
                else:
                    label_list.append("E = " + split_out[-1])
        else:
            if 'lin' in split_out[1]:
                if split_out[2] == '1111':
                    label_list.append("Linear 0.1 < E < " + split_out[-1])
                else:
                    label_list.append("Linear " + split_out[2] + " < E < " + split_out[-1])
            if 'exp' in split_out[1]:
                label_list.append("Exp. " + split_out[2] + " < E < " + split_out[-1])
            else:
                if split_out[-1] == '1111':
                    label_list.append("E = 0.1 Old Model")
                else:
                    label_list.append("E = " + split_out[-1] + " Old Model")
    return label_list


file_in = sys.argv[1]
sphere_r = 10
if file_in=='big_sphere':
    sphere_r = 20
shrink_r = str(float(sys.argv[2]))
xy_res = 141.7/512
z_res = 120/151
xy_bound = 141.7
z_bound = 120
nbeads = 4000
color_list = ['r', 'orange', 'y', 'chartreuse', 'g', 'deepskyblue']
color_list = ['g', 'r', 'k']

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

tag_list = ['3.0_const_1111',
            '3.0_const_100',
            '3.0_lin_1111_100']

x_values = np.linspace((xy_bound/2) + sphere_r, xy_bound, 1000)
y = xy_bound/2
z = z_bound/2

fig, ax = plt.subplots()
iter_count = 0

for save_tag in tag_list:

    if iter_count==100:
        old_mesh = Mesh()
        with XDMFFile("../ChangingModulus/MeshesXDMF/sphere.xdmf") as infile:
            infile.read(old_mesh)
        mvc = MeshValueCollection("size_t", old_mesh, 3)
        with XDMFFile("../ChangingModulus/MeshesXDMF/sphere_function.xdmf") as infile:
            infile.read(mvc, "geometrical")
        old_subdomains = cpp.mesh.MeshFunctionSizet(old_mesh, mvc)

        gmark, cmark = 0, 100
        dx = Measure('dx', domain=old_mesh, subdomain_data=old_subdomains, metadata={'quadrature_degree': 2})

        V = VectorFunctionSpace(old_mesh, "Lagrange", 1)
        du = TrialFunction(V)
        v = TestFunction(V)
        old_u = Function(V)
        hdf5_file = HDF5File(old_mesh.mpi_comm(),
                             "../ChangingModulus/saved_functions/sphere_0.70_" + save_tag + "_func.h5",
                             "r")
        hdf5_file.read(old_u, "/function")
        hdf5_file.close()
        disp = np.array([np.linalg.norm(old_u(x, y, z)) for x in x_values])
    else:
        hdf5_file = HDF5File(mesh.mpi_comm(),
                             "saved_functions/" + file_in + "_" + save_tag + "_func.h5",
                             "r")
        hdf5_file.read(u, "/function")
        hdf5_file.close()
        disp = np.array([np.linalg.norm(u(x, y, z)) for x in x_values])
    zord=1
    aa = 1
    ls = '-'
    # if iter_count == 1:
    #     zord = 1
    #     aa = .5
    #     ls='-'
    # else:
    #     zord = 100
    #     aa = .5
    #     ls='--'
    ax.plot(x_values - xy_bound/2, disp,
            LineWidth=2,
            c=color_list[iter_count],
            zorder=zord,
            alpha=aa,
            linestyle=ls)
    iter_count = iter_count + 1

ax.legend(tags2labels(tag_list))
ax.set_ylabel('Displacement (um)')
ax.set_xlabel('Distance from center of cell (um)')
ax.set_title("Continuous Displacement")
ax.set_aspect(15)
ax.set_ylim([0, 3.5])
plt.show()
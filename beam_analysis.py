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
            elif 'exp' in split_out[1]:
                label_list.append("Exp. " + split_out[2] + " < E < " + split_out[-1])
            else:
                if split_out[-1] == '1111':
                    label_list.append("E = 0.1 Old Model")
                else:
                    label_list.append("E = " + split_out[-1] + " Old Model")
    return label_list


# file_in = sys.argv[1]
file_in = "beam_hq"

color_list = ['r', 'orange', 'y', 'chartreuse', 'g', 'deepskyblue']
color_list = ['g', 'b', 'k', 'g', 'b']

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

tag_list = ['1_lin_100_0.01',
            '1_exp_100_0.01',
            '1_const_1']

x_values = np.linspace(0, 100, 2000)
y = 5
z = 5

fig, ax = plt.subplots()
iter_count = 0

for save_tag in tag_list:

    hdf5_file = HDF5File(mesh.mpi_comm(),
                         "saved_functions/" + file_in + "_1.0_" + save_tag + "_func.h5",
                         "r")
    hdf5_file.read(u, "/function")
    hdf5_file.close()
    disp = np.array([np.linalg.norm(u(x, y, z)) for x in x_values])

    zord=1
    aa = 1
    ls = '-'
    if iter_count>2:
        ls = "--"
    ax.plot(np.flip(x_values, axis=0),
            disp,
            LineWidth=2,
            c=color_list[iter_count],
            zorder=zord,
            alpha=aa,
            linestyle=ls)
    iter_count = iter_count + 1

ax.legend(['Linear 0.01 < E < 100', 'Exp. 0.01 < E < 100', 'E = 1'])
ax.set_ylabel('Displacement (um)')
ax.set_xlabel('Distance from moving boundary (um)')
ax.set_title("Continuous Displacement")
# ax.set_aspect(.5)
ax.set_ylim([0, 1.5])
plt.savefig("pics_temp/temp.png")
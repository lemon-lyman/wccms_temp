import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import sys
import os
import meshio
import fileinput
from scipy.spatial import Delaunay



def remap_mesh(mesh):
    remap_start = time.time()

    points = mesh.points
    old_tetra = mesh.cells['tetra']
    print("old points shape: ", points.shape)
    print("old tet shape: ", old_tetra.shape)
    referenced_idx = np.unique(old_tetra)
    referenced_points = points[referenced_idx]
    new_idx = np.arange(referenced_idx.shape[0])
    mapping = dict(zip(referenced_idx, new_idx))
    new_tetra = np.vectorize(mapping.get)(old_tetra)
    mesh.points = referenced_points
    mesh.cells = {'tetra':new_tetra}
    print("new points shape: ", referenced_points.shape)
    print("new tet shape: ", new_tetra.shape)

    print("remap time: ", time.time() - remap_start)
    return mesh

old_mesh = meshio.read("MeshesMSH/beam_hq.msh")
new_mesh = remap_mesh(old_mesh)
meshio.write("MeshesXDMF/beam_hq.xdmf", meshio.Mesh(points=new_mesh.points, cells={"tetra": new_mesh.cells["tetra"]}))
meshio.write("MeshesXDMF/beam_hq_function.xdmf", meshio.Mesh(points=new_mesh.points, cells={"tetra": new_mesh.cells["tetra"]},
                        cell_data={"tetra": {"geometrical": new_mesh.cell_data["tetra"]["gmsh:geometrical"]}}))

# read_path = 'skeleton_stuff/'
# write_path = 'skeleton_stuff/'
# files = []
# temp_list = ['long.msh']
# for file in temp_list:
#     # if 'geometrical' not in file:
#     if True:
#         # if len(file) == 9:
#         if True:
#             print(file)
#             files.append(file)
#             new_mesh = remap_mesh(meshio.read(read_path + file))
#             meshio.write(write_path + file[:-4] + ".xdmf", meshio.Mesh(points=new_mesh.points, cells={"tetra": new_mesh.cells["tetra"]}))
#             meshio.write(write_path + file[:-4] + "_function.xdmf", meshio.Mesh(points=new_mesh.points, cells={"tetra": new_mesh.cells["tetra"]},
#                                     cell_data={"tetra": {"geometrical": new_mesh.cell_data["tetra"]["gmsh:geometrical"]}}))
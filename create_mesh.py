import numpy as np
import meshio
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


a = 20
b = a
c = 20

offset_arr = [141.7/2, 141.7/2, 120/2]

n = 30
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def foo(x_):
    return x_*2

def get_coords(theta, psi):
    return np.array([a * np.sin(theta) * np.cos(psi) + offset_arr[0],
                     b * np.sin(theta) * np.sin(psi) + offset_arr[1],
                     c * np.cos(theta) + offset_arr[2]])

theta_vec = np.linspace(0, np.pi, n)
psi_vec = np.linspace(0, 2*np.pi, n)

angle_comb = np.array(np.meshgrid(theta_vec, psi_vec)).T.reshape(-1, 2)

ellipse_points = np.asarray([get_coords(ang[0], ang[1]) for ang in angle_comb])
triangulation = ConvexHull(ellipse_points)
# ax.scatter(ellipse_points[:, 0],
#            ellipse_points[:, 1],
#            ellipse_points[:, 2])
# [ax.plot([ellipse_points[tri[0], 0], ellipse_points[tri[1], 0]],
#         [ellipse_points[tri[0], 1], ellipse_points[tri[1], 1]],
#         [ellipse_points[tri[0], 2], ellipse_points[tri[1], 2]], c='k') for tri in triangulation.simplices]
# [ax.plot([ellipse_points[tri[1], 0], ellipse_points[tri[2], 0]],
#         [ellipse_points[tri[1], 1], ellipse_points[tri[2], 1]],
#         [ellipse_points[tri[1], 2], ellipse_points[tri[2], 2]], c='k') for tri in triangulation.simplices]
# [ax.plot([ellipse_points[tri[2], 0], ellipse_points[tri[0], 0]],
#         [ellipse_points[tri[2], 1], ellipse_points[tri[0], 1]],
#         [ellipse_points[tri[2], 2], ellipse_points[tri[0], 2]], c='k') for tri in triangulation.simplices]
# plt.show()
mesh = meshio.Mesh(ellipse_points,
                   {'triangle': triangulation.simplices})
meshio.write("MeshesMSH/big_sphere.msh", mesh)
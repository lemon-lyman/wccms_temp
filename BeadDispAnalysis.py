import numpy as np
import meshio
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix as dist_mat


def get_bins(distances, displacements):
    dx = 10
    bins = np.arange(10, 101, dx)
    y = []
    n = []
    for b in bins:
        inds = np.where(np.logical_and(np.greater_equal(distances, b - dx),
                                np.less_equal(distances, b)))[0]
        n.append(inds.shape[0])
        y.append(displacements[inds].mean())
    return bins-(dx/2), y, n


shrink_list = [.9, .85, .8, .75, .7]
color_list = ['r', 'orange', 'y', 'chartreuse', 'g', 'deepskyblue']


cells = ['sphere', 'sphere', 'sphere', 'sphere', 'sphere', 'sphere', 'sphere', 'sphere', 'sphere']
tags = ['const_0.01', 'const_0.000001', 'const_0.00001', 'const_0.0001', 'const_0.001', 'const_1111', 'const_1', 'const_10', 'const_100']
count = 0
for cell, tag in zip(cells, tags):
    cell_surface = meshio.read('MeshesMSH/' + cell + '.msh')
    cell_points = cell_surface.points
    fig, ax = plt.subplots()
    for idx, shrink in enumerate(shrink_list[::-1]):

        ## Normal loading
        # beads = np.loadtxt("beads_from_func/" + cell + "_" + str(shrink).ljust(4, '0') + "_" + tag + "_beads.txt", delimiter=",")
        # disps = np.loadtxt("beads_from_func/" + cell + "_" + str(shrink).ljust(4, '0') + "_" + tag + "_disps.txt", delimiter=",")
        ## Loading for check_randomness
        beads = np.loadtxt("simultaneous_beads/" + cell + "_beads.txt", delimiter=",")
        disps = np.loadtxt("simultaneous_beads/" + cell + "_" + str(shrink).ljust(4, '0') + "_" + tag + "_disps.txt", delimiter=",")

        disp_mag = np.linalg.norm(disps, axis=1)
        d_mat = dist_mat(beads, cell_points)
        dist2cell = d_mat.min(1)
        # fig, ax = plt.subplots()
        # ax.scatter(dist2cell, disp_mag, c='k', s=3)
        xx, yy, n = get_bins(dist2cell, disp_mag)
        ax.plot(xx, yy, '-o', c=color_list[idx])

        if idx==0:
            temp_x_lim = ax.get_xlim()
            temp_y_lim = ax.get_ylim()

            #########################
            ax2 = ax.twinx()
            ax2.bar(xx + idx - 2, n, 1, color = color_list[idx], zorder=-1, alpha=.1)
            ax2.set_ylabel('n')
            #########################

        else:

            #########################
            ax2.bar(xx + idx - 2, n, 1, color = color_list[idx], zorder=-1, alpha=.1)
            #########################

            ax.set_xlim(temp_x_lim)
            ax.set_ylim(temp_y_lim)
    if count==0:
        overall_x_lim = ax.get_xlim()
        overall_y_lim = ax.get_ylim()

    else:
        ax.set_xlim(overall_x_lim)
        ax.set_ylim(overall_y_lim)

    #########################
    ax2.set_ylim([0, 1000])
    #########################

    split_tag = tag.split('_')
    if '1111' in split_tag:
        split_tag[split_tag.index('1111')] = '10e-1'

    if len(split_tag)==3:
        if split_tag[0] == 'lin':
            ax.set_title(cell + " Linear " + split_tag[1] + " < E < " + split_tag[2])
        else:
            ax.set_title(cell + " Exponential " + split_tag[1] + " < E < " + split_tag[2])
    else:
        if split_tag[-1] == '0.01':
            ax.set_title(cell + " E = 10e-2")
        elif split_tag[-1] == '0.001':
            ax.set_title(cell + " E = 10e-3")
        elif split_tag[-1] == '0.0001':
            ax.set_title(cell + " E = 10e-4")
        elif split_tag[-1] == '0.00001':
            ax.set_title(cell + " E = 10e-5")
        elif split_tag[-1] == '0.000001':
            ax.set_title(cell + " E = 10e-6")
        elif split_tag[-1] == '1':
            ax.set_title(cell + " E = 10e0")
        elif split_tag[-1] == '10':
            ax.set_title(cell + " E = 10e1")
        elif split_tag[-1] == '100':
            ax.set_title(cell + " E = 10e2")
        else:
            ax.set_title(cell + " E = 10e-1")




    ax.set_xlabel("Distance From Cell (um)")
    ax.set_ylabel("Bead Displacement (um)")
    ax.legend([str(s).ljust(4, '0') for s in shrink_list[::-1]])
    plt.show()
    count = count + 1

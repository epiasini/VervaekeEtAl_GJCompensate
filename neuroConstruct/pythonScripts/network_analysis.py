import pdb
import csv
import numpy as np

from scipy.spatial import distance
from matplotlib import pyplot as plt

if __name__ == "__main__":
    n_trials = 10
    positions = np.loadtxt(open("cell_positions.csv", "rb"), delimiter=",")
    positions = positions.reshape(n_trials, positions.shape[0]/n_trials, 3)
    distances = np.array([distance.pdist(p) for p in positions])
    distances_squareform = np.array([distance.squareform(v) for v in distances])

    # load edge list for each trial
    edges = np.loadtxt(open("edge_lists.csv", "rb"), delimiter=",", dtype=np.int)
    edge_lists = []
    previous_edge = (None, None)
    for e in edges:
        if e[0] == 0 and previous_edge[0] != 0:
            edge_lists.append([])
        edge_lists[-1].append(e)
        previous_edge = e
    for k,l in enumerate(edge_lists):
        edge_lists[k] = np.array(l)
    # for every trial, create a temporary list where all edges are ordered (i<j)
    ordered_edge_lists = [l.tolist() for l in edge_lists]
    for edge_list in ordered_edge_lists:
        for k,e in enumerate(edge_list):
            if e[0] > e[1]:
                edge_list[k] = [e[1], e[0]]
    ordered_edge_lists = ((tuple(l2) for l2 in l1) for l1 in ordered_edge_lists)
    # for every trial, filter the ordered edge list and keep only its
    # unique elements. The result is a list of lists of all connected
    # cell pairs, regardless of the number of gap junctions between
    # them
    cell_pair_lol = [np.array(list(set(l))) for l in ordered_edge_lists]
    # for each trial, for each connected pair, calculate the distance
    # between somata
    edge_lengths = []
    for trial, cell_pair_list in enumerate(cell_pair_lol):
        lengths = np.array([distances_squareform[trial][tuple(e)] for e in cell_pair_list])
        edge_lengths.append(lengths)

    print distances[0].shape
    print edge_lists[0].shape
    print cell_pair_lol[0].shape
    print edge_lengths[0].shape

    hist_all = np.histogram(distances, range=[0.,160.])
    hist_connected = np.histogram(np.concatenate(edge_lengths), range=[0.,160.])

    fig, ax1 = plt.subplots()
    ax1.hist(x=(np.ravel(distances), np.concatenate(edge_lengths)), range=[0., 160.])
    ax2 = ax1.twinx()
    print hist_all[1]
    x = hist_all[1][:-1]+(hist_all[1][1]-hist_all[1][0])/2
    y = -1745 + 1836/(1 + np.exp((x-267)/39))
    y = y/y.max()
    ax2.plot(x, np.nan_to_num(np.asarray(hist_connected[0], dtype=np.float)/np.asarray(hist_all[0], dtype=np.float)), marker='o', color='r')
    ax2.plot(x, y, marker='o', color='k')
    ax2.set_ylim((0., 1.1))
    plt.show()

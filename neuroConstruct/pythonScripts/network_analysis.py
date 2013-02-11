import csv
import numpy as np

from scipy.spatial import distance
from scipy import optimize
from matplotlib import pyplot as plt

def vervaeke2010_spatial_dependence(r):
    return (-1745. + 1836./(1 + np.exp((r-267)/39)))/100

def degree_integrand(r, r_0, c):
    return r/(1+np.exp((r-r_0)/c))

def approximate_analytical_integral(r_max, r_0, c):
    return c*r_0*(r_max/c - np.log(1+np.exp(r_max/c)) + np.log(2)) + (c*np.pi)**2/12

def infinite_net_analytical_sol(rho_goc=4.6e-6, l=80., a=0.856, r_0=122, delta=16.9):
    return rho_goc * 2 * np.pi * l * a * (r_0**2/2 + (delta*np.pi)**2/12)

def calculate_degree(n=28, r_max=155., a=-17.45, b=-18.36, c=39, r_0=267, dx=0.01):
    x = np.arange(0, r_max, dx)
    numerical_sol = (n-1)*(a - (2*b)*np.trapz(degree_integrand(x, r_0=r_0, c=c), dx=dx)/r_max**2)
    print np.trapz(degree_integrand(x, r_0=r_0, c=c), dx=dx)
    print approximate_analytical_integral(r_max,r_0,c)
    approximate_analytical_sol = (n-1)*(a - (2*b)*approximate_analytical_integral(r_max,r_0,c)/r_max**2)
    return numerical_sol, infinite_net_analytical_sol()

def fermi_dirac_dependence(r, a, r_0, delta):
    return a/(1 + np.exp((r-r_0)/delta))

def fermi_dirac_fit():
    x = np.arange(0, 155., 0.01)
    vervaeke_values = vervaeke2010_spatial_dependence(x)
    a, r_0, delta = optimize.curve_fit(fermi_dirac_dependence, x, vervaeke_values, [0.8, 267, 39])[0]
    print a, r_0, delta
    new_x = np.arange(0, 200, 0.01)
    plt.plot(x, vervaeke_values)
    plt.plot(new_x, fermi_dirac_dependence(new_x, a, r_0, delta))
    plt.grid('on')
    plt.show()

if __name__ == "__main__":
    n_trials = 10
    positions = np.loadtxt(open("cell_positions.csv", "rb"), delimiter=",")
    n_cells = positions.shape[0]/n_trials
    positions = positions.reshape(n_trials, n_cells, 3)
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
    # for each trial, create a temporary list where all edges are ordered (i<j)
    ordered_edge_lists = [l.tolist() for l in edge_lists]
    for edge_list in ordered_edge_lists:
        for k,e in enumerate(edge_list):
            if e[0] > e[1]:
                edge_list[k] = [e[1], e[0]]
    ordered_edge_lists = ((tuple(l2) for l2 in l1) for l1 in ordered_edge_lists)
    # for each trial, filter the ordered edge list and keep only its
    # unique elements. The result is a list of lists of all connected
    # cell pairs, regardless of the number of gap junctions between
    # them
    cell_pair_lol = [np.array(list(set(l))) for l in ordered_edge_lists]
    # for each trial, calculate the degree sequence of the
    # network. Note that the degree of a cell is defined as the number
    # of unique cells it's connected to; the number of gap junctions
    # on a cell is an obvious upper bound for the cell's degree.
    degree_sequences = np.array([[len([pair for pair in l if any((pair[0]==k, pair[1]==k))]) for k in range(n_cells)] for l in cell_pair_lol])
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
    print degree_sequences
    print degree_sequences.mean()

    hist_all = np.histogram(distances, range=[0.,160.])
    hist_connected = np.histogram(np.concatenate(edge_lengths), range=[0.,160.])

    fig, ax1 = plt.subplots()
    ax1.hist(x=(np.ravel(distances), np.concatenate(edge_lengths)), range=[0., 160.])
    ax2 = ax1.twinx()
    # figure out the centres of the histogram bars on the x axis
    x = hist_all[1][:-1]+(hist_all[1][1]-hist_all[1][0])/2
    # plot the conected/total cells ratio
    ax2.plot(x, np.nan_to_num(np.asarray(hist_connected[0], dtype=np.float)/np.asarray(hist_all[0], dtype=np.float)), marker='o', color='r')
    # plot the fit to the experimental data (Vervaeke2010, figure 7)
    ax2.plot(x, vervaeke2010_spatial_dependence(x), marker='o', color='k')
    ax2.set_ylim((0., 1.1))
    plt.show()

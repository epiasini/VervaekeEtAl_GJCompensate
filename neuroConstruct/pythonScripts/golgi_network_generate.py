'''
Generate n_trials different Golgi network with the Vervaeke2012 algorithm, and save cell positions and edge lists to disk as csv files.
'''
import time
import random
import math
import csv

def euclidean_distance(u,v):
    return math.sqrt(sum([(u[k]-v[k])**2 for k in range(len(u))]))

if __name__ == "__main__":
    from java.lang import System, Float, Long
    from java.io import File
    from java.util import Vector, ArrayList
    from ucl.physiol import neuroconstruct as nc

    n_trials = 10

    pm = nc.project.ProjectManager(None,None)
    project_path = '../VervaekeEtAl-GJCompensate.ncx'
    project_file = File(project_path)
    project = pm.loadProject(project_file)

    sim_config_name = 'Default Simulation Configuration'
    sim_config = project.simConfigInfo.getSimConfig(sim_config_name)

    cell_positions_file = open('cell_positions.csv', 'wb')
    edge_lists_file = open('edge_lists.csv', 'wb')

    cell_positions_writer = csv.writer(cell_positions_file)
    edge_lists_writer = csv.writer(edge_lists_file)

    for trial in range(n_trials):
        nC_seed = Long(random.getrandbits(32))
        # generate
        pm.doGenerate(sim_config_name, nC_seed)
        while pm.isGenerating():
            time.sleep(0.02)
        print('network generated')

        n_cells = project.generatedCellPositions.getNumberInAllCellGroups()
        syn_conn = project.generatedNetworkConnections.getSynapticConnections('NetConn_CellGroup_4_CellGroup_4_1')

        cell_positions = [(pos_record.x_pos, pos_record.y_pos, pos_record.z_pos) for pos_record in project.generatedCellPositions.getPositionRecords('CellGroup_4')]
        edges = [(sc.sourceEndPoint.cellNumber, sc.targetEndPoint.cellNumber) for sc in syn_conn if sc.props[0].weight]
        #pair_distances = [euclidean_distance(cell_positions[i], cell_positions[j]) for i,j in edges]
        #degree_sequence = [float(len([e for e in edges if any((e[0]==n, e[1]==n))])) for n in range(n_cells)]
        #mean_degree = sum(degree_sequence)/len(degree_sequence)

        cell_positions_writer.writerows(cell_positions)
        edge_lists_writer.writerows(edges)




    System.exit(0)


def test_euclidean_distance():
    import numpy as np
    for n_dims in range(1,10):
        u = np.random.random(size=n_dims) - 0.5
        v = np.random.random(size=n_dims) - 0.5
        assert euclidean_distance(u,v) == np.linalg.norm(u - v)



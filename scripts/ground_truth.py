#! /usr/bin/env python
# encoding: utf-8

import numpy as np
import sys
import trimesh
import networkx as nx
import os


mesh = trimesh.load_mesh(sys.argv[1])
mesh_name = os.path.splitext(os.path.basename(sys.argv[1]))[0]

graph = trimesh.graph.vertex_adjacency_graph(mesh)

geodesics = nx.shortest_path_length(graph, source=0, weight='weight')

N = len(mesh.vertices)

ground_truth = np.zeros((N))
period = 2*np.pi / 0.3 * 2

for i in range(N):
    ground_truth[i] = np.sin(0.2*geodesics.get(i))

np.savetxt("rsc/truth/" + mesh_name + "_truth.csv", ground_truth)
np.savetxt("rsc/truth/" + mesh_name + "_vertices.csv", mesh.vertices)
np.savetxt("rsc/truth/" + mesh_name + "_faces.csv", mesh.faces)

# mesh.show()

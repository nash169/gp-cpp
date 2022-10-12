#! /usr/bin/env python
# encoding: utf-8

import numpy as np
import sys
import trimesh
import networkx as nx
import os

name = sys.argv[1] if len(sys.argv) > 1 else "sphere"

# Load mesh
mesh = trimesh.load_mesh(
    trimesh.interfaces.gmsh.load_gmsh("rsc/" + name + ".msh"))

# edges without duplication
edges = mesh.edges_unique

# the actual length of each unique edge
length = mesh.edges_unique_length

# create the graph with edge attributes for length
graph = nx.Graph()
for edge, L in zip(edges, length):
    graph.add_edge(*edge, length=L)

# graph = trimesh.graph.vertex_adjacency_graph(mesh)

geodesics = nx.shortest_path_length(graph, source=0, weight='length')

N = len(mesh.vertices)

ground_truth = np.zeros((N))
period = 2  # 2*np.pi / 0.3 * 2

for i in range(N):
    ground_truth[i] = 2 * np.sin(geodesics.get(i) * period + 0.3)

if not os.path.exists("outputs"):
    os.mkdir("outputs")

if not os.path.exists("outputs/truth"):
    os.mkdir("outputs/truth")

np.savetxt("outputs/truth/" + name + "_truth.csv", ground_truth)
np.savetxt("outputs/truth/" + name + "_vertices.csv", mesh.vertices)
np.savetxt("outputs/truth/" + name + "_faces.csv", mesh.faces)

# mesh.show()

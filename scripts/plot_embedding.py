#! /usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import sys

from matplotlib import cm
from io_utils import get_data

# Data
mesh_name = sys.argv[1] if len(sys.argv) > 1 else "sphere"
num_modes = int(sys.argv[2]) if len(sys.argv) > 2 else 5

# Nodes
nodes = np.loadtxt("outputs/truth/"+mesh_name+"_vertices.csv")

# Samples
samples = np.loadtxt("outputs/truth/"+mesh_name+"_reference.csv")

# Fem modes
data = get_data("outputs/modes/fem_"+mesh_name+"_eigs.000000", "eigs")
eigs_fem = data["eigs"][:num_modes]

modes_fem = {}
data = get_data("outputs/modes/fem_"+mesh_name+"_modes.000000", "modes")
for i in range(len(eigs_fem)):
    modes_fem[i] = data["modes"][:, i]

# Diffusion modes
data = get_data("outputs/modes/diffusion_"+mesh_name+"_eigs.000000", "eigs")
eigs_diffusion = data["eigs"][:num_modes]

modes_diffusion = {}
data = get_data("outputs/modes/diffusion_"+mesh_name+"_modes.000000", "modes")
for i in range(len(eigs_diffusion)):
    modes_diffusion[i] = data["modes"][:, i]

# Free data
del data

# Plot nodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((np.ptp(nodes[:, 0]), np.ptp(
    nodes[:, 1]), np.ptp(nodes[:, 2])))
ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2])
ax.set_title('Sampled nodes from manifold')
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$x_3$')

# Plot nodes + reference
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((np.ptp(nodes[:, 0]), np.ptp(
    nodes[:, 1]), np.ptp(nodes[:, 2])))
ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2])
ax.scatter(samples[:, 0], samples[:, 1], samples[:, 2], s=20, c="red", alpha=1)
ax.set_title('Sampled nodes as ground truth')
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$x_3$')

# Plot fem spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(range(len(eigs_fem)), eigs_fem, '-o')
ax.axis('square')
ax.set_title('FEM Spectrum')
ax.set_xlabel('$i$')
ax.set_ylabel('$\lambda$')

# Plot diffusion spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(range(len(eigs_diffusion)), eigs_diffusion, '-o')
ax.axis('square')
ax.set_title('Diffusion Maps Spectrum')
ax.set_xlabel('$i$')
ax.set_ylabel('$\lambda$')

# Plot fem embedding
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((np.ptp(modes_fem[1]), np.ptp(
    modes_fem[2]), np.ptp(modes_fem[3])))
ax.scatter(modes_fem[1], modes_fem[2], modes_fem[3])
ax.set_title('FEM Embedding')
ax.set_xlabel('$\psi_1$')
ax.set_ylabel('$\psi_2$')
ax.set_zlabel('$\psi_3$')

# Plot diffusion embedding
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((np.ptp(modes_diffusion[1]), np.ptp(
    modes_diffusion[2]), np.ptp(modes_diffusion[3])))
ax.scatter(modes_diffusion[1], modes_diffusion[2], modes_diffusion[3])
ax.set_title('Diffusion Maps Embedding')
ax.set_xlabel('$\psi_1$')
ax.set_ylabel('$\psi_2$')
ax.set_zlabel('$\psi_3$')

plt.show()

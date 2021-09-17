#! /usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import sys

from matplotlib import cm
from io_utils import get_data

data = get_data(sys.argv[1], "SAMPLES", "FUNCTION", "EIGVEC")

res = 100

fig = plt.figure()
ax = fig.add_subplot(121, projection="3d")
surf = ax.plot_surface(
    data["SAMPLES"][:, 0].reshape((res, res), order='F'),
    data["SAMPLES"][:, 1].reshape((res, res), order='F'),
    data["SAMPLES"][:, 2].reshape((res, res), order='F'),
    facecolors=cm.jet(data["FUNCTION"].reshape(
        (res, res), order='F')/np.amax(data["FUNCTION"])),
    antialiased=True, linewidth=0
)
fig.colorbar(surf, ax=ax)

# x.quiver(data["X"], data["Y"], data["F"].reshape(data["X"].shape, order='F'),
#           data["G"][:, 0].reshape(data["X"].shape, order='F'), data["G"][:, 1].reshape(
#               data["X"].shape, order='F'), data["G"][:, 2].reshape(data["X"].shape, order='F'),
#           length=0.1, color='black')


ax = fig.add_subplot(122)
contour = ax.contourf(
    data["SAMPLES"][:, 0].reshape((res, res), order='F'),
    data["SAMPLES"][:, 1].reshape((res, res), order='F'),
    data["SAMPLES"][:, 2].reshape((res, res), order='F'),
    cmap=cm.jet,
    antialiased=True, linewidth=0
)
ax.set_aspect("equal", "box")
fig.colorbar(contour, ax=ax)

# ax.quiver(data["X"], data["Y"],
#           data["G"][:, 0].reshape(data["X"].shape, order='F'), data["G"][:, 1].reshape(
#               data["X"].shape, order='F'))

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
# ax.set_box_aspect((np.ptp(data["SAMPLES"][:, 0]), np.ptp(
#     data["SAMPLES"][:, 1]), np.ptp(data["SAMPLES"][:, 2])))
ax.scatter(data["SAMPLES"][:, 0], data["SAMPLES"][:, 1], data["SAMPLES"][:, 2])
ax.set_title('Sampled nodes from manifold')
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$x_3$')

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.set_box_aspect((np.ptp(data["EIGVEC"][:, 1]), np.ptp(
    data["EIGVEC"][:, 2]), np.ptp(data["EIGVEC"][:, 3])))
ax.scatter(data["EIGVEC"][:, 1], data["EIGVEC"][:, 2], data["EIGVEC"][:, 3])
ax.set_title('Embedding')
ax.set_xlabel('$u_1$')
ax.set_ylabel('$u_2$')
ax.set_zlabel('$u_3$')

plt.show()

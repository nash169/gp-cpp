#! /usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import sys

from matplotlib import cm
from io_utils import get_data

# data = get_data(sys.argv[1], "X", "Y", "sample_x", "sample_y", "GP", "OPT")

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(data["X"], data["Y"], color='black', linestyle='dashed')
# ax.scatter(data["sample_x"], data["sample_y"], color="red")
# ax.plot(data["X"], data["GP"], color="green")
# ax.plot(data["X"], data["OPT"], color="blue")

# Number of extracted modes
num_modes = 5

# Get modes
modes_fem = {}
modes_diffusion = {}
for i in range(num_modes):
    data = get_data("rsc/modes/fem_sphere_mode_" +
                    str(i)+".000000", "Ordering: 0")
    modes_fem[i] = data["Ordering: 0"]
    data = get_data("rsc/modes/diffusion_sphere_mode_" +
                    str(i)+".000000", "mode")
    modes_diffusion[i] = data["mode"]

fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
ax.set_box_aspect((np.ptp(modes_fem[1]), np.ptp(
    modes_fem[2]), np.ptp(modes_fem[3])))
ax.scatter(modes_fem[1], modes_fem[2], modes_fem[3])

ax = fig.add_subplot(122, projection='3d')
ax.set_box_aspect((np.ptp(modes_diffusion[1]), np.ptp(
    modes_diffusion[2]), np.ptp(modes_diffusion[3])))
ax.scatter(modes_diffusion[1], modes_diffusion[2], modes_diffusion[3])

plt.show()

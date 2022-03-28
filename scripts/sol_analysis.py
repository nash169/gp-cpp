#! /usr/bin/env python
# encoding: utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
from io_utils import get_data

mesh_name = sys.argv[1] if len(sys.argv) > 1 else "sphere"
alg_name = sys.argv[2] if len(sys.argv) > 2 else "ambient"

mse = get_data("rsc/solutions/"+alg_name+"_"+mesh_name+"_gp.csv", "mse")["mse"]

mean = np.mean(mse, axis=0)
std = np.std(mse, axis=0)
iter_pace = np.array([25, 50, 75, 100, 125, 150])

colors = np.array(["#377eb8", "#ff7f00", "#4daf4a", "#f781bf",
                   "#a65628", "#984ea3", "#999999", "#e41a1c", "#dede00", ])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.fill_between(iter_pace,
                mean - std,
                mean + std,
                alpha=0.25,
                linewidth=0,
                color=colors[0])
ax.plot(iter_pace, mean, linewidth=2, color=colors[1])
ax.set_xlabel("# Training points", fontsize=15)
ax.set_ylabel("MSE", fontsize=15)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.tick_params(axis="x", direction="out")
ax.tick_params(axis="y", length=0)
for spine in ax.spines.values():
    spine.set_position(("outward", 5))
ax.set_axisbelow(True)
ax.grid(axis="y", color="0.9", linestyle="-", linewidth=1)

fig.savefig("rsc/solutions/"+alg_name+"_"+mesh_name+"_analysis.png")
# plt.show()

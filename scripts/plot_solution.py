#! /usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import sys

from matplotlib import cm
from io_utils import get_data

data = get_data("outputs/solutions/1d_gp.csv", "X", "Y",
                "sample_x", "sample_y", "GP", "OPT", "VARIANCE")

colors = ["#377eb8", "#ff7f00", "#4daf4a", "#f781bf",
          "#a65628", "#984ea3", "#999999", "#e41a1c", "#dede00"]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(data["X"], data["Y"], color='black', linestyle='dashed')
ax.scatter(data["sample_x"], data["sample_y"], color="red")
ax.plot(data["X"], data["GP"], color="green")
ax.plot(data["X"], data["OPT"], color="blue")
ax.fill_between(data["X"], data["OPT"] - 3*data["VARIANCE"], data["OPT"] + 3*data["VARIANCE"],
                alpha=0.25, linewidth=0, color=colors[1])

plt.show()

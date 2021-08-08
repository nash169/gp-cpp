#! /usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import sys

from matplotlib import cm
from io_utils import get_data

data = get_data("rsc/solutions/1d_gp.csv", "X", "Y",
                "sample_x", "sample_y", "GP", "OPT")

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(data["X"], data["Y"], color='black', linestyle='dashed')
ax.scatter(data["sample_x"], data["sample_y"], color="red")
ax.plot(data["X"], data["GP"], color="green")
ax.plot(data["X"], data["OPT"], color="blue")

plt.show()

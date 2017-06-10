#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np

from mpl_toolkits.mplot3d import Axes3D


def draw_3D_histogramm(color_list, output_file):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y = np.random.rand(2, 100) * 4
    hist, xedges, yedges = np.histogram2d(x, y, bins=4)

    elements = (len(xedges) - 1) * (len(yedges) - 1)
    xpos, ypos = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.25)

    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)
    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = hist.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')

    plt.savefig( "tratata.svg", facecolor="#D6D6D6")
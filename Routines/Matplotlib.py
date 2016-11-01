__author__ = 'mahajrod'
import os
import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector, BboxConnectorPatch
from matplotlib.lines import Line2D

import numpy as np


class MatplotlibRoutines:
    def __init__(self):
        pass

    @staticmethod
    def connect_bbox(bbox1, bbox2,
                     loc1a, loc2a, loc1b, loc2b,
                     prop_lines, prop_patches=None):
        if prop_patches is None:
            prop_patches = prop_lines.copy()
            prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

        c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
        c1.set_clip_on(False)
        c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
        c2.set_clip_on(False)

        bbox_patch1 = BboxPatch(bbox1, **prop_patches)
        bbox_patch2 = BboxPatch(bbox2, **prop_patches)

        p = BboxConnectorPatch(bbox1, bbox2,
                               loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                               **prop_patches)
        p.set_clip_on(False)

        return c1, c2, bbox_patch1, bbox_patch2, p

    def zoom_effect(self, ax1, ax2, xmin, xmax, alpha=0.22, color="gray", **kwargs):
        """
        ax1 : the main axes
        ax2 : the zoomed axes
        (xmin,xmax) : the limits of the colored area in both plot axes.

        connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
        be marked.  The keywords parameters will be used ti create
        patches.

        """

        trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
        trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

        bbox = Bbox.from_extents(xmin, 0, xmax, 1)

        mybbox1 = TransformedBbox(bbox, trans1)
        mybbox2 = TransformedBbox(bbox, trans2)

        prop_patches = kwargs.copy()
        prop_patches["ec"] = "none"
        prop_patches["alpha"] = alpha
        prop_patches["color"] = color

        c1, c2, bbox_patch1, bbox_patch2, p = self.connect_bbox(mybbox1, mybbox2,
                                                                loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                                                                prop_lines=kwargs, prop_patches=prop_patches)

        ax1.add_patch(bbox_patch1)
        ax2.add_patch(bbox_patch2)
        ax2.add_patch(c1)
        ax2.add_patch(c2)
        ax2.add_patch(p)

        return c1, c2, bbox_patch1, bbox_patch2, p

    @staticmethod
    def percent_histogram(data, output_prefix, n_bins=20, title="", xlabel="%", ylabel="Number", label=None,
                          extensions=("png", "svg"), legend=None, legend_location="best", input_mode="percent", xmax=None,
                          xmin=None):

        figure = plt.figure()
        subplot = plt.subplot(1, 1, 1)

        plt.hist(data, bins=n_bins, label=label)
        if input_mode == "percent":
            plt.xlim(xmin=0, xmax=100)
        elif input_mode == "fraction":
            plt.xlim(xmin=0, xmax=1)
        else:
            raise ValueError("Unrecognized type of input data(neither percents nor fractions)")
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if legend:
            plt.legend((legend,), loc=legend_location)
        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

    def percent_histogram_from_file(self, data_file, output_prefix, data_type=float, column_list=None, separator=None,
                                    comments="#", n_bins=20, title="", xlabel="%", ylabel="Number",
                                    extensions=("png", "svg"), legend=None, legend_location="best",
                                    stats_as_legend=False, input_mode="percent", label=None):
        #print column_list
        data = np.loadtxt(data_file, dtype=data_type, comments=comments, delimiter=separator, usecols=column_list)
        if input_mode == "percent":
            n_bins = np.linspace(0, 100, n_bins+1)
        elif input_mode == "fraction":
            n_bins = np.linspace(0, 1.0, n_bins+1)
        else:
            raise ValueError("Unrecognized type of input data(neither percents nor fractions)")
        legenda = "Total: %i\nMean: %.2f %%\nMedian: %.2f %%" if input_mode == "percent" else "Total: %i\nMean: %.2f\nMedian: %.2f"
        legenda = legenda % (len(data), np.mean(data), np.median(data)) if stats_as_legend else legend
        
        self.percent_histogram(data, output_prefix=output_prefix, n_bins=n_bins, title=title, xlabel=xlabel,
                               ylabel=ylabel, extensions=extensions, legend=legenda, legend_location=legend_location,
                               input_mode=input_mode, label=label)

    @staticmethod
    def add_line(axes, start, end, color):
        line = Line2D([start[0], end[0]], [start[1], end[1]], color=color)
        return axes.add_line(line)

    @staticmethod
    def int_histogram(data, output_prefix, n_bins=None, title="", xlabel="%", ylabel="Number",
                      extensions=("png", "svg"), legend=None, legend_location="best"):



        figure = plt.figure()
        subplot = plt.subplot(1, 1, 1)

        maximum = np.max(data)
        minimum = np.min(data)
        number_of_bins = n_bins if n_bins else maximum - minimum + 1
        plt.hist(data, bins=number_of_bins)

        plt.xlim(xmin=minimum - 1, xmax=maximum)

        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if legend:
            plt.legend((legend,), loc=legend_location)
        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

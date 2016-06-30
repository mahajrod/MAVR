#!/usr/bin/env python
__author__ = 'mahajrod'

from collections import OrderedDict
from matplotlib.patches import Rectangle


class FeatureStyle(object):
    def __init__(self, edgecolor=None, facecolor=None, main_color=None, visible=True):
        self.edgecolor = edgecolor
        self.facecolor = facecolor
        self.main_color = main_color
        self.visible = visible


class RectangularFeatureStyle(FeatureStyle):

    def __init__(self, xy, width, height, angle=0.0, subfeatures=OrderedDict(), **kwargs):
        """
        angle: rotation in degrees(anti-clockwise)
        """
        super(RectangularFeatureStyle, self).__init__(xy, width, height, angle=angle, **kwargs)
        self.subfeatures = subfeatures



class ExonStyle():
    def __init__(self, height, color, border_color, border_width, strandness=False):
        self.height = height
        self.color = color
        self.border_color = border_color
        self.border_width = border_width
        self.strandness = strandness


class IntronStyle():
    def __init__(self):
        pass


class CDSStyle():
    def __init__(self):
        pass


class UTRStyle():
    def __init__(self):
        pass


class GeneStyle():
    def __init__(self, intron_style, exon_style, cds_style=None, five_utr_style=None, three_utr_style=None):
        self.intron = intron_style
        self.exon = exon_style
        self.five_utr = five_utr_style
        self.cds = cds_style
        self.three_utr = three_utr_style


class DomenStyle():
    def __init__():
        pass
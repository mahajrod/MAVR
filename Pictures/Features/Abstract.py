
from collections import OrderedDict
from matplotlib.patches import Rectangle


class Feature(object):
    def __init__(self, edgecolor=None, facecolor=None, main_color=None, visible=True):
        self.edgecolor = edgecolor
        self.facecolor = facecolor
        self.main_color = main_color
        self.visible = visible


class RectangularFeature(Rectangle, object):
    """
    Inherits from matplotlib.patches.Rectangle(which is based on matplotlib.patches.Patch)
    """
    def __init__(self, xy, width, height, angle=0.0, subfeatures=OrderedDict(), feature_type=None,
                 visible_features_list=[], **kwargs):
        """
        angle: rotation in degrees(anti-clockwise)
        """
        super(RectangularFeature, self).__init__(xy, width, height, angle=angle, **kwargs)
        self.subfeatures = subfeatures
        self.feature_type = feature_type
        self.visible_features_list = visible_features_list


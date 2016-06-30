__author__ = 'mahajrod'
from collections import OrderedDict

from Pictures.Features.Abstract import RectangularFeature


class RectangularProtein(RectangularFeature):
    def __init__(self, xy, width, height, angle=0.0, subfeatures=OrderedDict(), visible_features_list=[], **kwargs):
        """
        angle: rotation in degrees(anti-clockwise)
        """
        super(RectangularProtein, self).__init__(xy, width, height, subfeatures=subfeatures,
                                                 feature_type="protein", angle=angle,
                                                 visible_features_list=visible_features_list if visible_features_list else ["domen"], **kwargs)
#!/usr/bin/env python

from Pipelines.Abstract import Pipeline


class GapClosingPipeline(Pipeline):
    def __init__(self):
        Pipeline.__init__(self)

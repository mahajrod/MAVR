#!/usr/bin/env python

from Pipelines.Abstract import Pipeline


class ContigAssemblyPipeline(Pipeline):
    def __init__(self):
        Pipeline.__init__(self)

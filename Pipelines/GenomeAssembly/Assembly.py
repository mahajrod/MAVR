#!/usr/bin/env python

from Pipelines.Abstract import Pipeline
from Pipelines.GenomeAssembly.GapClosing import GapClosingPipeline
from Pipelines.GenomeAssembly.Scaffolding import ScaffoldingPipeline
from Pipelines.GenomeAssembly.ContigAssembly import ContigAssemblyPipeline


class AssemblyPipeline(Pipeline, GapClosingPipeline, ScaffoldingPipeline, ContigAssemblyPipeline):
    def __init__(self):
        Pipeline.__init__(self)

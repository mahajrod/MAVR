#!/usr/bin/env python

from Pipelines.GenomeAssembly.ContigAssembly import ContigAssemblyPipeline
from Pipelines.GenomeAssembly.Scaffolding import ScaffoldingPipeline
from Pipelines.GenomeAssembly.GapClosing import GapClosingPipeline
from Pipelines.GenomeAssembly.Assembly import AssemblyPipeline

ContigAssemblyPipeline = ContigAssemblyPipeline()
ScaffoldingPipeline = ScaffoldingPipeline()
GapClosingPipeline = GapClosingPipeline()
AssemblyPipeline = AssemblyPipeline()

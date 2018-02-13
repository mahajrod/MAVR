#!/usr/bin/env python
#import os
#import shutil
#from collections import OrderedDict
#from Tools.Filter import Cookiecutter, Trimmomatic, FaCut
from Tools.Abstract import Tool
from Routines.Matplotlib import MatplotlibRoutines
from Routines.Fastq import FastQRoutines
#from Parsers.FaCut import FaCutReport
#from Parsers.Coockiecutter import CoockiecutterReport
#from Parsers.Trimmomatic import TrimmomaticReport

#from CustomCollections.GeneralCollections import TwoLvlDict


class Pipeline(Tool, MatplotlibRoutines, FastQRoutines):
    def __init__(self, max_threads=1, max_memory=10):
        Tool.__init__(self, cmd="", max_threads=max_threads, max_memory=max_memory)

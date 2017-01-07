#!/usr/bin/env python
#import os
#import shutil
#from collections import OrderedDict
#from Tools.Filter import Cookiecutter, Trimmomatic, FaCut
from Tools.Abstract import Tool
from Routines.Matplotlib import MatplotlibRoutines
#from Parsers.FaCut import FaCutReport
#from Parsers.Coockiecutter import CoockiecutterReport
#from Parsers.Trimmomatic import TrimmomaticReport

#from CustomCollections.GeneralCollections import TwoLvlDict


class Pipeline(Tool, MatplotlibRoutines):
    def __init__(self):
        Tool.__init__(self, cmd="")
